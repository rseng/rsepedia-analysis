![Logo](https://raw.githubusercontent.com/biosimulations/biosimulations/dev/libs/shared/assets/src/assets/images/biosimulations-logo/logo-white.svg)

[![App Status](https://deployment.biosimulations.org/api/badge?name=biosimulations-dev&revision=true)](https://deployment.biosimulations.org/applications/biosimulations-dev)
[![App Status](https://deployment.biosimulations.org/api/badge?name=biosimulations-prod&revision=true)](https://deployment.biosimulations.org/applications/biosimulations-prod)
[![Continuous Integration](https://github.com/biosimulations/biosimulations/workflows/Continuous%20Integration/badge.svg)](https://github.com/biosimulations/biosimulations/actions?query=workflow%3A%22Continuous+Integration%22)
[![Continuous Deployment](https://github.com/biosimulations/biosimulations/workflows/Continuous%20Deployment/badge.svg)](https://github.com/biosimulations/biosimulations/actions?query=workflow%3A%22Continuous+Deployment%22)

[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/5204/badge)](https://bestpractices.coreinfrastructure.org/projects/5204)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Maintainability Rating](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=sqale_rating)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Reliability Rating](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=reliability_rating)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Security Rating](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=security_rating)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)

[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=bugs)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Code Smells](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=code_smells)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Duplicated Lines (%)](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=duplicated_lines_density)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=ncloc)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Vulnerabilities](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=vulnerabilities)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)
[![Technical Debt](https://sonarcloud.io/api/project_badges/measure?project=biosimulations_biosimulations&metric=sqale_index)](https://sonarcloud.io/summary/new_code?id=biosimulations_biosimulations)



[![All Contributors](https://img.shields.io/github/all-contributors/biosimulations/biosimulations/HEAD)](#contributors-)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-1.4-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![Commitizen friendly](https://img.shields.io/badge/commitizen-friendly-brightgreen.svg)](http://commitizen.github.io/cz-cli/)

[![DOI](https://zenodo.org/badge/207730765.svg)](https://zenodo.org/badge/latestdoi/207730765) 

# BioSimulations 🧬

More comprehensive and more predictive models have the potential to advance biology, bioengineering, and medicine. Building more predictive models will likely require the collaborative efforts of many investigators. This requires teams to be able to share and reuse model components and simulations. Despite extensive efforts to develop standards such as [COMBINE/OMEX](https://combinearchive.org/), [SBML](http://sbml.org), and [SED-ML](https://sed-ml.org), it remains difficult to reuse many models and simulations. One challenge to reusing models and simulations is the diverse array of incompatible modeling formats and simulation tools.

This package provides three tools which address this challenge:

- [BioSimulators](https://biosimulators.org) is a registry of containerized simulation tools that provide consistent interfaces. BioSimulators makes it easier to find and run simulations.
- [runBioSimulations](https://run.biosimulations.org) is a simple web application for using the BioSimulators containers to run simulations. This tool makes it easy to run a broad range of simulations without having to install any software.
- [BioSimulations](https://biosimulations.org) is a platform for sharing and running modeling studies. BioSimulations provides a central place for investigators to exchange studies. BioSimulations uses the BioSimulators simulation tools, and builds on the functionality of runBioSimulations.

This package provides the code for the BioSimulations, runBioSimulations, and BioSimulations websites, as well as the code for the backend services for all three applications. The package is implemented in TypeScript using Angular, NestJS, MongoDB, and Mongoose.

## Getting started ▶️

### Users 💻

Please use the hosted versions of BioSimulations, runBioSimulations, and BioSimulators at [https://biosimulations.org](https://biosimulations.org), [https://run.biosimulations.org](https://run.biosimulations.org), and [https://biosimulators.org](https://biosimulators.org).

Tutorials, help and information can be found at [https://docs.biosimulations.org](https://docs.biosimulations.org)

### Developers 🖥️

We welcome contributions to BioSimulations, runBioSimulations, and BioSimulations! Please see the [developer guide](https://docs.biosimulations.org/developers) for information about how to get started including how to install this package and how to run BioSimulations, runBioSimulations, and BioSimulators locally.

## License ⚖️

This package is released under the [MIT license](./License.md). This package uses a number of open-source third-party packages. Their licenses are summarized in [Dependencies](./about/Dependencies).

## Show your support 🤝

If you find this project interesting or useful, please give our repo a ⭐ and share with others that may benefit. If you use the code and tools in this repository as a part of an academic work, please cite us using the following bibtex entry. 

```
@software{Shaikh_BioSimulations,
author = {Shaikh, Bilal and Marupilla, Gnaneswara and Wilson, Mike and Michael, Blinov L. and Moraru, Ion I. and Karr, Jonathan R.},
doi = {10.5281/zenodo.5057108},
license = {MIT},
title = {{BioSimulations}},
url = {https://github.com/biosimulations/biosimulations}
}
```

## Contributors 🧑‍🤝‍🧑

This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York and the [Center for Cell Analysis and Modeling](https://health.uconn.edu/cell-analysis-modeling/) at UConn Health as part of the [Center for Reproducible Biomodeling Modeling](https://reproduciblebiomodels.org).

Numerous individuals and groups have contributed to BioSimulations, including:     

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/AMICI-dev"><img src="https://avatars.githubusercontent.com/u/68919097?v=4?s=100" width="100px;" alt=""/><br /><sub><b>AMICI</b></sub></a><br /><a href="#tool-AMICI-dev" title="Tools">🔧</a></td>
    <td align="center"><a href="https://fun.bio.keio.ac.jp/"><img src="https://avatars.githubusercontent.com/u/1589676?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Akira Funahashi</b></sub></a><br /><a href="#tool-funasoul" title="Tools">🔧</a></td>
    <td align="center"><a href="https://hellix.com/Alan/"><img src="https://avatars.githubusercontent.com/u/602265?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Alan Garny</b></sub></a><br /><a href="#ideas-agarny" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-agarny" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/ajelenak"><img src="https://avatars.githubusercontent.com/u/7267124?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Aleksandar Jelenak</b></sub></a><br /><a href="#tool-ajelenak" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/ASinanSaglam"><img src="https://avatars.githubusercontent.com/u/11724447?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ali Sinan Saglam</b></sub></a><br /><a href="#data-ASinanSaglam" title="Data">🔣</a></td>
    <td align="center"><a href="https://uni-tuebingen.de/en/127116"><img src="https://avatars.githubusercontent.com/u/1740827?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Andreas Dräger</b></sub></a><br /><a href="#tool-draeger" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/AnkitaxPriya"><img src="https://avatars.githubusercontent.com/u/44089458?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ankita</b></sub></a><br /><a href="#data-AnkitaxPriya" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://ankursinha.in/"><img src="https://avatars.githubusercontent.com/u/102575?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ankur Sinha</b></sub></a><br /><a href="#tool-sanjayankur31" title="Tools">🔧</a></td>
    <td align="center"><a href="https://research.pasteur.fr/en/member/anna-zhukova"><img src="https://avatars.githubusercontent.com/u/10465838?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Anna Zhukova</b></sub></a><br /><a href="#data-annazhukova" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/AnneGoelzer"><img src="https://avatars.githubusercontent.com/u/32333634?v=4?s=100" width="100px;" alt=""/><br /><sub><b>AnneGoelzer</b></sub></a><br /><a href="#data-AnneGoelzer" title="Data">🔣</a></td>
    <td align="center"><a href="https://www.mountsinai.org/profiles/arthur-p-goldberg"><img src="https://avatars.githubusercontent.com/u/33882?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Arthur P Goldberg</b></sub></a><br /><a href="#ideas-artgoldberg" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="http://aurelien.naldi.info/"><img src="https://avatars.githubusercontent.com/u/250984?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Aurélien Naldi</b></sub></a><br /><a href="#data-aurelien-naldi" title="Data">🔣</a> <a href="#tool-aurelien-naldi" title="Tools">🔧</a></td>
    <td align="center"><a href="http://bshaikh.com"><img src="https://avatars.githubusercontent.com/u/32490144?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Bilal Shaikh</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=bilalshaikh42" title="Code">💻</a> <a href="https://github.com/biosimulations/biosimulations/commits?author=bilalshaikh42" title="Documentation">📖</a> <a href="#infra-bilalshaikh42" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
    <td align="center"><a href="https://www.ebi.ac.uk/biomodels"><img src="https://avatars.githubusercontent.com/u/74367888?v=4?s=100" width="100px;" alt=""/><br /><sub><b>BioModels</b></sub></a><br /><a href="#data-EBI-BioModels" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://teusinkbruggemanlab.nl/brett-olivier/"><img src="https://avatars.githubusercontent.com/u/5011985?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brett Olivier</b></sub></a><br /><a href="#tool-bgoli" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/briandrawert"><img src="https://avatars.githubusercontent.com/u/1413538?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brian Drawert</b></sub></a><br /><a href="#tool-briandrawert" title="Tools">🔧</a></td>
    <td align="center"><a href="https://briansimulator.org/"><img src="https://avatars.githubusercontent.com/u/2292949?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brian simulator</b></sub></a><br /><a href="#tool-brian-team" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.copasi.org/"><img src="https://avatars.githubusercontent.com/u/1854399?v=4?s=100" width="100px;" alt=""/><br /><sub><b>COPASI</b></sub></a><br /><a href="#tool-copasi" title="Tools">🔧</a></td>
    <td align="center"><a href="https://reproduciblebiomodels.org"><img src="https://avatars.githubusercontent.com/u/70044163?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Center for Reproducible Biomedical Modeling</b></sub></a><br /><a href="#financial-reproducible-biomedical-modeling" title="Financial">💵</a> <a href="#fundingFinding-reproducible-biomedical-modeling" title="Funding Finding">🔍</a> <a href="#projectManagement-reproducible-biomedical-modeling" title="Project Management">📆</a></td>
    <td align="center"><a href="https://github.com/CiaranWelsh"><img src="https://avatars.githubusercontent.com/u/19502680?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ciaran Welsh</b></sub></a><br /><a href="#tool-CiaranWelsh" title="Tools">🔧</a></td>
    <td align="center"><a href="https://claudine-chaouiya.pedaweb.univ-amu.fr/index.html"><img src="https://avatars.githubusercontent.com/u/40125033?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Claudine Chaouiya</b></sub></a><br /><a href="#data-chaouiya" title="Data">🔣</a> <a href="#tool-chaouiya" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/danv61"><img src="https://avatars.githubusercontent.com/u/29076329?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dan Vasilescu</b></sub></a><br /><a href="#tool-danv61" title="Tools">🔧</a></td>
    <td align="center"><a href="https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/5122/index.html"><img src="https://avatars.githubusercontent.com/u/18048784?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Daniel Weindl</b></sub></a><br /><a href="#tool-dweindl" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/dbrnz"><img src="https://avatars.githubusercontent.com/u/239220?v=4?s=100" width="100px;" alt=""/><br /><sub><b>David Brooks</b></sub></a><br /><a href="#tool-dbrnz" title="Tools">🔧</a></td>
    <td align="center"><a href="http://about.me/david.nickerson"><img src="https://avatars.githubusercontent.com/u/811244?v=4?s=100" width="100px;" alt=""/><br /><sub><b>David Nickerson</b></sub></a><br /><a href="#ideas-nickerso" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/DeepaMahm"><img src="https://avatars.githubusercontent.com/u/29662579?v=4?s=100" width="100px;" alt=""/><br /><sub><b>DeepaMahm</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/issues?q=author%3ADeepaMahm" title="Bug reports">🐛</a></td>
    <td align="center"><a href="https://github.com/jdieg0"><img src="https://avatars.githubusercontent.com/u/6570972?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Diego</b></sub></a><br /><a href="#tool-jdieg0" title="Tools">🔧</a></td>
    <td align="center"><a href="https://dilawars.me/"><img src="https://avatars.githubusercontent.com/u/895681?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dilawar Singh</b></sub></a><br /><a href="#tool-dilawar" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/edkerk"><img src="https://avatars.githubusercontent.com/u/7326655?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Eduard Kerkhoven</b></sub></a><br /><a href="#tool-edkerk" title="Tools">🔧</a></td>
    <td align="center"><a href="https://eagmon.github.io/"><img src="https://avatars.githubusercontent.com/u/6809431?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Eran Agmon</b></sub></a><br /><a href="#ideas-eagmon" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/Ermentrout"><img src="https://avatars.githubusercontent.com/u/7952422?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ermentrout</b></sub></a><br /><a href="#tool-Ermentrout" title="Tools">🔧</a></td>
    <td align="center"><a href="https://escher.github.io/"><img src="https://avatars.githubusercontent.com/u/9327950?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Escher</b></sub></a><br /><a href="#tool-escher" title="Tools">🔧</a></td>
    <td align="center"><a href="https://scholar.harvard.edu/fabianfroehlich/home"><img src="https://avatars.githubusercontent.com/u/14923969?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Fabian Fröhlich</b></sub></a><br /><a href="#tool-FFroehlich" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/zhangfengkai"><img src="https://avatars.githubusercontent.com/u/38113699?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Fengkai Zhang</b></sub></a><br /><a href="#tool-zhangfengkai" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/fbergmann"><img src="https://avatars.githubusercontent.com/u/949059?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Frank Bergmann</b></sub></a><br /><a href="#ideas-fbergmann" title="Ideas, Planning, & Feedback">🤔</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/GINsim"><img src="https://avatars.githubusercontent.com/u/32065286?v=4?s=100" width="100px;" alt=""/><br /><sub><b>GINsim</b></sub></a><br /><a href="#tool-GINsim" title="Tools">🔧</a></td>
    <td align="center"><a href="http://gmarupilla.com"><img src="https://avatars.githubusercontent.com/u/53095348?v=4?s=100" width="100px;" alt=""/><br /><sub><b>GMarupilla</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=gmarupilla" title="Code">💻</a></td>
    <td align="center"><a href="http://helikarlab.org/"><img src="https://avatars.githubusercontent.com/u/17307008?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Helikar Lab Personal</b></sub></a><br /><a href="#tool-HelikarLabPersonal" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.sys-bio.org/"><img src="https://avatars.githubusercontent.com/u/1054990?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Herbert Sauro</b></sub></a><br /><a href="#ideas-hsauro" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/hsorby"><img src="https://avatars.githubusercontent.com/u/778048?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Hugh Sorby</b></sub></a><br /><a href="#tool-hsorby" title="Tools">🔧</a></td>
    <td align="center"><a href="http://identifiers.org/"><img src="https://avatars.githubusercontent.com/u/18701545?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Idenfitiers.org</b></sub></a><br /><a href="#data-identifiers-org" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/JanHasenauer"><img src="https://avatars.githubusercontent.com/u/12297214?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jan Hasenauer</b></sub></a><br /><a href="#tool-JanHasenauer" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://bionetgen.org/"><img src="https://avatars.githubusercontent.com/u/8277248?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jim Faeder</b></sub></a><br /><a href="#tool-jrfaeder" title="Tools">🔧</a> <a href="#data-jrfaeder" title="Data">🔣</a></td>
    <td align="center"><a href="http://vcell.org"><img src="https://avatars.githubusercontent.com/u/20616724?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jim Schaff</b></sub></a><br /><a href="#ideas-jcschaff" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/jmrohwer"><img src="https://avatars.githubusercontent.com/u/502289?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Johann Rohwer</b></sub></a><br /><a href="#tool-jmrohwer" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/jhgennari"><img src="https://avatars.githubusercontent.com/u/2684850?v=4?s=100" width="100px;" alt=""/><br /><sub><b>John Gennari</b></sub></a><br /><a href="#ideas-jhgennari" title="Ideas, Planning, & Feedback">🤔</a> <a href="#tool-jhgennari" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/jreadey"><img src="https://avatars.githubusercontent.com/u/7785492?v=4?s=100" width="100px;" alt=""/><br /><sub><b>John Readey</b></sub></a><br /><a href="#tool-jreadey" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/johnsekar"><img src="https://avatars.githubusercontent.com/u/1610689?v=4?s=100" width="100px;" alt=""/><br /><sub><b>John Sekar</b></sub></a><br /><a href="#ideas-johnsekar" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/joncison"><img src="https://avatars.githubusercontent.com/u/1506863?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jon Ison</b></sub></a><br /><a href="#data-joncison" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://www.karrlab.org"><img src="https://avatars.githubusercontent.com/u/2848297?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jonathan Karr</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=jonrkarr" title="Code">💻</a> <a href="https://github.com/biosimulations/biosimulations/commits?author=jonrkarr" title="Documentation">📖</a> <a href="#design-jonrkarr" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/jtcooper10"><img src="https://avatars.githubusercontent.com/u/42880781?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Joshua Cooper</b></sub></a><br /><a href="#tool-jtcooper10" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/starboerg"><img src="https://avatars.githubusercontent.com/u/5522086?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jörn Starruß</b></sub></a><br /><a href="#tool-starboerg" title="Tools">🔧</a></td>
    <td align="center"><a href="http://juergen.pahle.de/"><img src="https://avatars.githubusercontent.com/u/5473011?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jürgen Pahle</b></sub></a><br /><a href="#tool-jpahle" title="Tools">🔧</a></td>
    <td align="center"><a href="https://www.karrlab.org/"><img src="https://avatars.githubusercontent.com/u/13785824?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Karr whole-cell modeling lab</b></sub></a><br /><a href="#ideas-KarrLab" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/0u812"><img src="https://avatars.githubusercontent.com/u/7402146?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Kyle Medley</b></sub></a><br /><a href="#tool-0u812" title="Tools">🔧</a> <a href="#ideas-0u812" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="http://lems.github.io/LEMS"><img src="https://avatars.githubusercontent.com/u/3033237?v=4?s=100" width="100px;" alt=""/><br /><sub><b>LEMS</b></sub></a><br /><a href="#tool-LEMS" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://loicpauleve.name/"><img src="https://avatars.githubusercontent.com/u/228657?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Loïc Paulevé</b></sub></a><br /><a href="#data-pauleve" title="Data">🔣</a> <a href="#tool-pauleve" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/luciansmith"><img src="https://avatars.githubusercontent.com/u/1736150?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Lucian Smith</b></sub></a><br /><a href="#ideas-luciansmith" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/lutzbrusch"><img src="https://avatars.githubusercontent.com/u/13622401?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Lutz Brusch</b></sub></a><br /><a href="#tool-lutzbrusch" title="Tools">🔧</a></td>
    <td align="center"><a href="https://uk.linkedin.com/in/manuelbernal"><img src="https://avatars.githubusercontent.com/u/8855107?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Manuel Bernal Llinares</b></sub></a><br /><a href="#data-mbdebian" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/MarcDinh"><img src="https://avatars.githubusercontent.com/u/50445930?v=4?s=100" width="100px;" alt=""/><br /><sub><b>MarcDinh</b></sub></a><br /><a href="#tool-MarcDinh" title="Tools">🔧</a></td>
    <td align="center"><a href="https://livermetabolism.com/"><img src="https://avatars.githubusercontent.com/u/900538?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Matthias König</b></sub></a><br /><a href="#ideas-matthiaskoenig" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://orcid.org/0000-0002-1509-4981"><img src="https://avatars.githubusercontent.com/u/992660?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Matúš Kalaš</b></sub></a><br /><a href="#data-matuskalas" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/vcellmike"><img src="https://avatars.githubusercontent.com/u/29076280?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michael Blinov</b></sub></a><br /><a href="#ideas-vcellmike" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="http://dumontierlab.com/"><img src="https://avatars.githubusercontent.com/u/993852?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michel Dumontier</b></sub></a><br /><a href="#data-micheldumontier" title="Data">🔣</a></td>
    <td align="center"><a href="http://www.cds.caltech.edu/~mhucka"><img src="https://avatars.githubusercontent.com/u/1450019?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mike Hucka</b></sub></a><br /><a href="#tool-mhucka" title="Tools">🔧</a></td>
    <td align="center"><a href="https://hpc.uchc.edu"><img src="https://avatars.githubusercontent.com/u/400595?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mike Wilson</b></sub></a><br /><a href="#infra-mpw6" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
    <td align="center"><a href="http://modeldb.yale.edu/"><img src="https://avatars.githubusercontent.com/u/38667483?v=4?s=100" width="100px;" alt=""/><br /><sub><b>ModelDB</b></sub></a><br /><a href="#data-ModelDBRepository" title="Data">🔣</a></td>
    <td align="center"><a href="https://unseenbio.com/"><img src="https://avatars.githubusercontent.com/u/135653?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Moritz E. Beber</b></sub></a><br /><a href="#tool-Midnighter" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.nibib.nih.gov/"><img src="https://avatars.githubusercontent.com/u/12418167?v=4?s=100" width="100px;" alt=""/><br /><sub><b>National Institute of Biomedical Imaging and Bioengineering</b></sub></a><br /><a href="#financial-NIBIB" title="Financial">💵</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://nih.gov/"><img src="https://avatars.githubusercontent.com/u/52710462?v=4?s=100" width="100px;" alt=""/><br /><sub><b>National Institutes of Health</b></sub></a><br /><a href="#financial-NIHGOV" title="Financial">💵</a></td>
    <td align="center"><a href="https://nsf.gov/"><img src="https://avatars.githubusercontent.com/u/23663503?v=4?s=100" width="100px;" alt=""/><br /><sub><b>National Science Foundation</b></sub></a><br /><a href="#financial-NSF-open" title="Financial">💵</a></td>
    <td align="center"><a href="https://docs.neuroml.org/"><img src="https://avatars.githubusercontent.com/u/2727519?v=4?s=100" width="100px;" alt=""/><br /><sub><b>NeuroML</b></sub></a><br /><a href="#tool-NeuroML" title="Tools">🔧</a></td>
    <td align="center"><a href="http://neurosimlab.org/"><img src="https://avatars.githubusercontent.com/u/14202113?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Neurosim lab</b></sub></a><br /><a href="#tool-Neurosim-lab" title="Tools">🔧</a></td>
    <td align="center"><a href="https://opencor.ws/"><img src="https://avatars.githubusercontent.com/u/754570?v=4?s=100" width="100px;" alt=""/><br /><sub><b>OpenCOR</b></sub></a><br /><a href="#tool-opencor" title="Tools">🔧</a></td>
    <td align="center"><a href="https://models.physiomeproject.org/"><img src="https://avatars.githubusercontent.com/u/1114929?v=4?s=100" width="100px;" alt=""/><br /><sub><b>PMR2 - the software behind the Auckland Physiome Repository</b></sub></a><br /><a href="#data-PMR2" title="Data">🔣</a></td>
    <td align="center"><a href="http://www.opensourcebrain.org/"><img src="https://avatars.githubusercontent.com/u/1556687?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Padraig Gleeson</b></sub></a><br /><a href="#data-pgleeson" title="Data">🔣</a> <a href="#tool-pgleeson" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/Paytonco"><img src="https://avatars.githubusercontent.com/u/7064808?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Payton Thomas</b></sub></a><br /><a href="#tool-Paytonco" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.comp-sys-bio.org/"><img src="https://avatars.githubusercontent.com/u/2159130?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Pedro Mendes</b></sub></a><br /><a href="#tool-pmendes" title="Tools">🔧</a></td>
    <td align="center"><a href="http://pedromonteiro.org/"><img src="https://avatars.githubusercontent.com/u/2027375?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Pedro T. Monteiro</b></sub></a><br /><a href="#tool-ptgm" title="Tools">🔧</a></td>
    <td align="center"><a href="http://pysces.sourceforge.net/"><img src="https://avatars.githubusercontent.com/u/6103247?v=4?s=100" width="100px;" alt=""/><br /><sub><b>PySCeS: The Python Simulator for Cellular Systems, provides a variety of tools for the analysis of cellular systems</b></sub></a><br /><a href="#tool-PySCeS" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/Ragzz1995"><img src="https://avatars.githubusercontent.com/u/16513966?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Raghul Kannan</b></sub></a><br /><a href="#tool-Ragzz1995" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/rsmsheriff"><img src="https://avatars.githubusercontent.com/u/7849690?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Rahuman Sheriff</b></sub></a><br /><a href="#data-rsmsheriff" title="Data">🔣</a></td>
    <td align="center"><a href="https://raashika03.github.io/rashika.rathi/"><img src="https://avatars.githubusercontent.com/u/45493793?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Rashika Rathi</b></sub></a><br /><a href="#data-raashika03" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://allencell.org/"><img src="https://avatars.githubusercontent.com/u/9079?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ryan Spangler</b></sub></a><br /><a href="#ideas-prismofeverything" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/Ryannjordan"><img src="https://avatars.githubusercontent.com/u/86376602?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ryann Jordan</b></sub></a><br /><a href="#data-Ryannjordan" title="Data">🔣</a></td>
    <td align="center"><a href="http://sbml.org/About"><img src="https://avatars.githubusercontent.com/u/1799692?v=4?s=100" width="100px;" alt=""/><br /><sub><b>SBML Team</b></sub></a><br /><a href="#tool-sbmlteam" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/SED-ML"><img src="https://avatars.githubusercontent.com/u/29736746?v=4?s=100" width="100px;" alt=""/><br /><sub><b>SED-ML</b></sub></a><br /><a href="#tool-SED-ML" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/skeating"><img src="https://avatars.githubusercontent.com/u/1736558?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Sarah Keating</b></sub></a><br /><a href="#tool-skeating" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/shoops"><img src="https://avatars.githubusercontent.com/u/1760522?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Stefan Hoops</b></sub></a><br /><a href="#tool-shoops" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.smoldyn.org/"><img src="https://avatars.githubusercontent.com/u/33039297?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Steve Andrews</b></sub></a><br /><a href="#data-ssandrews" title="Data">🔣</a> <a href="#tool-ssandrews" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/StochSS"><img src="https://avatars.githubusercontent.com/u/3344600?v=4?s=100" width="100px;" alt=""/><br /><sub><b>StochSS</b></sub></a><br /><a href="#tool-StochSS" title="Tools">🔧</a></td>
    <td align="center"><a href="https://maiage.inrae.fr/en/biosys"><img src="https://avatars.githubusercontent.com/u/32363627?v=4?s=100" width="100px;" alt=""/><br /><sub><b>SysBioINRAe</b></sub></a><br /><a href="#tool-SysBioInra" title="Tools">🔧</a> <a href="#data-SysBioInra" title="Data">🔣</a></td>
    <td align="center"><a href="https://science.vu.nl/en/research/molecular-cell-biology/systems-bioinformatics/index.aspx"><img src="https://avatars.githubusercontent.com/u/12168054?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Systems Biology Lab, Vrije Universiteit Amsterdam</b></sub></a><br /><a href="#tool-SystemsBioinformatics" title="Tools">🔧</a></td>
    <td align="center"><a href="http://systemsbiology.ucsd.edu/"><img src="https://avatars.githubusercontent.com/u/4237829?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Systems Biology Research Group</b></sub></a><br /><a href="#data-SBRG" title="Data">🔣</a></td>
    <td align="center"><a href="https://www.hdfgroup.org/"><img src="https://avatars.githubusercontent.com/u/8572050?v=4?s=100" width="100px;" alt=""/><br /><sub><b>The HDF Group</b></sub></a><br /><a href="#tool-HDFGroup" title="Tools">🔧</a></td>
    <td align="center"><a href="http://neuron.yale.edu/"><img src="https://avatars.githubusercontent.com/u/38567601?v=4?s=100" width="100px;" alt=""/><br /><sub><b>The NEURON Simulator</b></sub></a><br /><a href="#tool-neuronsimulator" title="Tools">🔧</a></td>
    <td align="center"><a href="https://cellml.org/"><img src="https://avatars.githubusercontent.com/u/2141414?v=4?s=100" width="100px;" alt=""/><br /><sub><b>The home of CellML on Github</b></sub></a><br /><a href="#tool-cellml" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/metatoaster"><img src="https://avatars.githubusercontent.com/u/372914?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Tommy Yu</b></sub></a><br /><a href="#data-metatoaster" title="Data">🔣</a></td>
    <td align="center"><a href="https://www.itersdesktop.com/"><img src="https://avatars.githubusercontent.com/u/663341?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Tung Nguyen</b></sub></a><br /><a href="#data-ntung" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/sys-bio"><img src="https://avatars.githubusercontent.com/u/5590646?v=4?s=100" width="100px;" alt=""/><br /><sub><b>UW Sauro Lab</b></sub></a><br /><a href="#tool-sys-bio" title="Tools">🔧</a></td>
    <td align="center"><a href="https://vega.github.io/"><img src="https://avatars.githubusercontent.com/u/11796929?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Vega</b></sub></a><br /><a href="#tool-vega" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/veitveit"><img src="https://avatars.githubusercontent.com/u/15800709?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Veit Schwämmle</b></sub></a><br /><a href="#data-veitveit" title="Data">🔣</a></td>
    <td align="center"><a href="http://vcell.org/"><img src="https://avatars.githubusercontent.com/u/29076025?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Virtual Cell</b></sub></a><br /><a href="#tool-virtualcell" title="Tools">🔧</a> <a href="#data-virtualcell" title="Data">🔣</a></td>
    <td align="center"><a href="http://genome.jouy.inra.fr/~wliebermeis/index_en.html"><img src="https://avatars.githubusercontent.com/u/3976679?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Wolfram Liebermeister</b></sub></a><br /><a href="#tool-liebermeister" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/YinHoon"><img src="https://avatars.githubusercontent.com/u/11270172?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Yin Hoon Chew</b></sub></a><br /><a href="#ideas-YinHoon" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://www.linkedin.com/in/zakandrewking/"><img src="https://avatars.githubusercontent.com/u/1250400?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Zachary A. King</b></sub></a><br /><a href="#tool-zakandrewking" title="Tools">🔧</a> <a href="#data-zakandrewking" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/abulovic"><img src="https://avatars.githubusercontent.com/u/1510530?v=4?s=100" width="100px;" alt=""/><br /><sub><b>abulovic</b></sub></a><br /><a href="#data-abulovic" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/cjmyers"><img src="https://avatars.githubusercontent.com/u/3507191?v=4?s=100" width="100px;" alt=""/><br /><sub><b>cjmyers</b></sub></a><br /><a href="#ideas-cjmyers" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/dczielinski"><img src="https://avatars.githubusercontent.com/u/4442307?v=4?s=100" width="100px;" alt=""/><br /><sub><b>dczielinski</b></sub></a><br /><a href="#tool-dczielinski" title="Tools">🔧</a></td>
    <td align="center"><a href="https://thesustainablevegan.org/"><img src="https://avatars.githubusercontent.com/u/60083977?v=4?s=100" width="100px;" alt=""/><br /><sub><b>freiburgermsu</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=freiburgermsu" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/jtyurkovich"><img src="https://avatars.githubusercontent.com/u/5396263?v=4?s=100" width="100px;" alt=""/><br /><sub><b>jtyurkovich</b></sub></a><br /><a href="#tool-jtyurkovich" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://fun.bio.keio.ac.jp/software/libsbmlsim/"><img src="https://avatars.githubusercontent.com/u/16151392?v=4?s=100" width="100px;" alt=""/><br /><sub><b>libsbmlsim</b></sub></a><br /><a href="#tool-libsbmlsim" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/moraru"><img src="https://avatars.githubusercontent.com/u/7397814?v=4?s=100" width="100px;" alt=""/><br /><sub><b>moraru</b></sub></a><br /><a href="#infra-moraru" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#ideas-moraru" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/obodeit"><img src="https://avatars.githubusercontent.com/u/38722594?v=4?s=100" width="100px;" alt=""/><br /><sub><b>obodeit</b></sub></a><br /><a href="#tool-obodeit" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/opencobra"><img src="https://avatars.githubusercontent.com/u/2708410?v=4?s=100" width="100px;" alt=""/><br /><sub><b>openCOBRA</b></sub></a><br /><a href="#tool-opencobra" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/RuleWorld"><img src="https://avatars.githubusercontent.com/u/11491841?v=4?s=100" width="100px;" alt=""/><br /><sub><b>ruleworld</b></sub></a><br /><a href="#tool-RuleWorld" title="Tools">🔧</a> <a href="#data-RuleWorld" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/yexilein"><img src="https://avatars.githubusercontent.com/u/30040612?v=4?s=100" width="100px;" alt=""/><br /><sub><b>yexilein</b></sub></a><br /><a href="#tool-yexilein" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/z-haiman"><img src="https://avatars.githubusercontent.com/u/29131681?v=4?s=100" width="100px;" alt=""/><br /><sub><b>z-haiman</b></sub></a><br /><a href="#tool-z-haiman" title="Tools">🔧</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
A key to the above emojis is available [here](https://allcontributors.org/docs/en/emoji-key).


## Contributing to BioSimulations 🛠️
We enthusiastically welcome contributions to BioSimulations! Please see the [guide to contributing](docs/CONTRIBUTING.md) and the [developer's code of conduct](docs/CODE_OF_CONDUCT.md).

## Funding 💰

This package was developed with support from the National Institute for Bioimaging and Bioengineering (award P41EB023912).

## Questions and comments ❓

We welcome any comments, questions, or discussion about the project. Please create a discussion or question in our [discussion forum](https://github.com/biosimulations/biosimulations/discussions).

To privately contact the BioSimulations team, you can send us an email at [info@biosimulations.org](mailto:info@biosimulations.org).
# config-common

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test config-common` to execute the unit tests via [Jest](https://jestjs.io).
# config-angular

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test config-angular` to execute the unit tests.
# config-nest

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test config-nest` to execute the unit tests via [Jest](https://jestjs.io).
# datamodel-common

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test datamodel-common` to execute the unit tests via [Jest](https://jestjs.io).
# datamodel-database

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test datamodel-database` to execute the unit tests via [Jest](https://jestjs.io).
# datamodel-utils

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test datamodel-utils` to execute the unit tests via [Jest](https://jestjs.io).
# datamodel-api

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test datamodel-api` to execute the unit tests via [Jest](https://jestjs.io).
This folder should contain the API for the "core" objects of the api, corresponding to (but not necessarily exactly) the objects exposed by endpoints such as SimulationRuns and Specifications. 
# datamodel-simulation-runs

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test datamodel-simulation-runs` to execute the unit tests via [Jest](https://jestjs.io).
# auth-common

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test auth-common` to execute the unit tests via [Jest](https://jestjs.io).
# auth-angular

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test auth-angular` to execute the unit tests.
# Auth-Client

Connects to the auth0 service to get a M2M token to authenticate for calling APIs
The authentication information is read from the config service.
This client is only used for Machine to Machine authentication, and not for user accounts.
# auth-nest

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test auth-nest` to execute the unit tests via [Jest](https://jestjs.io).
# ontology-datamodel

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test ontology-datamodel` to execute the unit tests via [Jest](https://jestjs.io).
# ontology-utils

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test ontology-utils` to execute the unit tests via [Jest](https://jestjs.io).
# ontology-extra-sources

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test ontology-extra-sources` to execute the unit tests via [Jest](https://jestjs.io).
# ontology-sources

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test ontology-sources` to execute the unit tests via [Jest](https://jestjs.io).
# ontology-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test ontology-client` to execute the unit tests.
# ontology-api

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test ontology-api` to execute the unit tests via [Jest](https://jestjs.io).
# simulation-project-utils-service

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test simulation-project-utils-service` to execute the unit tests.
# simulation-project-utils-ui

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test simulation-project-utils-ui` to execute the unit tests.
# combine-api-angular-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test combine-api-angular-client` to execute the unit tests.
## @

### Building

To install the required dependencies and to build the typescript sources run:
```
npm install
npm run build
```

### publishing

First build the package then run ```npm publish dist``` (don't forget to specify the `dist` folder!)

### consuming

Navigate to the folder of your consuming project and run one of next commands.

_published:_

```
npm install @ --save
```

_without publishing (not recommended):_

```
npm install PATH_TO_GENERATED_PACKAGE/dist.tgz --save
```

_It's important to take the tgz file, otherwise you'll get trouble with links on windows_

_using `npm link`:_

In PATH_TO_GENERATED_PACKAGE/dist:
```
npm link
```

In your project:
```
npm link 
```

__Note for Windows users:__ The Angular CLI has troubles to use linked npm packages.
Please refer to this issue https://github.com/angular/angular-cli/issues/8284 for a solution / workaround.
Published packages are not effected by this issue.


#### General usage

In your Angular project:


```
// without configuring providers
import { ApiModule } from '';
import { HttpClientModule } from '@angular/common/http';

@NgModule({
    imports: [
        ApiModule,
        // make sure to import the HttpClientModule in the AppModule only,
        // see https://github.com/angular/angular/issues/20575
        HttpClientModule
    ],
    declarations: [ AppComponent ],
    providers: [],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```

```
// configuring providers
import { ApiModule, Configuration, ConfigurationParameters } from '';

export function apiConfigFactory (): Configuration {
  const params: ConfigurationParameters = {
    // set configuration parameters here.
  }
  return new Configuration(params);
}

@NgModule({
    imports: [ ApiModule.forRoot(apiConfigFactory) ],
    declarations: [ AppComponent ],
    providers: [],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```

```
// configuring providers with an authentication service that manages your access tokens
import { ApiModule, Configuration } from '';

@NgModule({
    imports: [ ApiModule ],
    declarations: [ AppComponent ],
    providers: [
      {
        provide: Configuration,
        useFactory: (authService: AuthService) => new Configuration(
          {
            basePath: environment.apiUrl,
            accessToken: authService.getAccessToken.bind(authService)
          }
        ),
        deps: [AuthService],
        multi: false
      }
    ],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```

```
import { DefaultApi } from '';

export class AppComponent {
    constructor(private apiGateway: DefaultApi) { }
}
```

Note: The ApiModule is restricted to being instantiated once app wide.
This is to ensure that all services are treated as singletons.

#### Using multiple OpenAPI files / APIs / ApiModules
In order to use multiple `ApiModules` generated from different OpenAPI files,
you can create an alias name when importing the modules
in order to avoid naming conflicts:
```
import { ApiModule } from 'my-api-path';
import { ApiModule as OtherApiModule } from 'my-other-api-path';
import { HttpClientModule } from '@angular/common/http';

@NgModule({
  imports: [
    ApiModule,
    OtherApiModule,
    // make sure to import the HttpClientModule in the AppModule only,
    // see https://github.com/angular/angular/issues/20575
    HttpClientModule
  ]
})
export class AppModule {

}
```


### Set service base path
If different than the generated base path, during app bootstrap, you can provide the base path to your service. 

```
import { BASE_PATH } from '';

bootstrap(AppComponent, [
    { provide: BASE_PATH, useValue: 'https://your-web-service.com' },
]);
```
or

```
import { BASE_PATH } from '';

@NgModule({
    imports: [],
    declarations: [ AppComponent ],
    providers: [ provide: BASE_PATH, useValue: 'https://your-web-service.com' ],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```


#### Using @angular/cli
First extend your `src/environments/*.ts` files by adding the corresponding base path:

```
export const environment = {
  production: false,
  API_BASE_PATH: 'http://127.0.0.1:8080'
};
```

In the src/app/app.module.ts:
```
import { BASE_PATH } from '';
import { environment } from '../environments/environment';

@NgModule({
  declarations: [
    AppComponent
  ],
  imports: [ ],
  providers: [{ provide: BASE_PATH, useValue: environment.API_BASE_PATH }],
  bootstrap: [ AppComponent ]
})
export class AppModule { }
```  
# combine-api-nest-client-wrapper

This library was generated with [Nx](https://nx.dev).


## Running unit tests

Run `nx test combine-api-nest-client-wrapper` to execute the unit tests via [Jest](https://jestjs.io).


# combine-api-nest-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test combine-api-nest-client` to execute the unit tests via [Jest](https://jestjs.io).
## @

### Building

To install the required dependencies and to build the typescript sources run:
```
npm install
npm run build
```

#### General usage

In your Nestjs project:


```
// without configuring providers
import { ApiModule } from '';
import { HttpModule } from '@nestjs/common';

@Module({
    imports: [
        ApiModule,
        HttpModule
    ],
    providers: []
})
export class AppModule {}
```

```
// configuring providers
import { ApiModule, Configuration, ConfigurationParameters } from '';

export function apiConfigFactory (): Configuration => {
  const params: ConfigurationParameters = {
    // set configuration parameters here.
  }
  return new Configuration(params);
}

@Module({
    imports: [ ApiModule.forRoot(apiConfigFactory) ],
    declarations: [ AppComponent ],
    providers: [],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```

```
import { DefaultApi } from '';

export class AppComponent {
    constructor(private apiGateway: DefaultApi) { }
}
```

Note: The ApiModule a dynamic module and instantiated once app wide.
This is to ensure that all services are treated as singletons.

#### Using multiple swagger files / APIs / ApiModules
In order to use multiple `ApiModules` generated from different swagger files,
you can create an alias name when importing the modules
in order to avoid naming conflicts:
```
import { ApiModule } from 'my-api-path';
import { ApiModule as OtherApiModule } from 'my-other-api-path';
import { HttpModule } from '@nestjs/common';

@Module({
  imports: [
    ApiModule,
    OtherApiModule,
    HttpModule
  ]
})
export class AppModule {

}
```


### Set service base path
If different than the generated base path, during app bootstrap, you can provide the base path to your service. 

```
import { BASE_PATH } from '';

bootstrap(AppComponent, [
    { provide: BASE_PATH, useValue: 'https://your-web-service.com' },
]);
```
or

```
import { BASE_PATH } from '';

@Module({
    imports: [],
    declarations: [ AppComponent ],
    providers: [ provide: BASE_PATH, useValue: 'https://your-web-service.com' ],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```


#### Using @nestjs/cli
First extend your `src/environments/*.ts` files by adding the corresponding base path:

```
export const environment = {
  production: false,
  API_BASE_PATH: 'http://127.0.0.1:8080'
};
```

In the src/app/app.module.ts:
```
import { BASE_PATH } from '';
import { environment } from '../environments/environment';

@Module({
  declarations: [
    AppComponent
  ],
  imports: [ ],
  providers: [
    { 
      provide: 'BASE_PATH', 
      useValue: environment.API_BASE_PATH 
    }
  ]
})
export class AppModule { }
```
# simulators-database-models

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test simulators-database-models` to execute the unit tests via [Jest](https://jestjs.io).
# shared-storage

A wrapper module for the NestJS Configures the connection to the S3 bucket used across the apps.

## Configuration

The configuration is loaded via the config module, and uses the following environment variables

- STORAGE_ENDPOINT - A URL for an AWS S3 compatible s3 serverr
- STORAGE_ACCESS - The Access Key used for connecting to the bucket
- STORAGE_SECRET - The Secret key associated with the key
- STORAGE_BUCKET - The name of the storage bucket to use

## Limitations

The module is marked as global to limit the initial setup that is needed. This contains a limitation of only being able to use one service account across a given app. If multiple accounts are needed, the S3 module can be configured independently for the app.
See [the S3 module documentation](https://github.com/svtslv/nestjs-s3#readme) for setup instructions
# shared-icons

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-icons` to execute the unit tests.
# shared-exceptions--exceptions

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test shared-exceptions--exceptions` to execute the unit tests via [Jest](https://jestjs.io).
# Exceptions

Exception handling for HTTP applications

## Filters

This library contains exception filters for NestJS applications. Each exception filter is added to the providers array of the [SharedExceptionsModule](./src/lib/shared-exceptions.module.ts). Importing the module to the App Module includes the filters. The filters set the Response body to [JSON-API Error Object](https://jsonapi.org/format/#error-objects)

## Exception Classes

The library also contains exception classes that are based on the [JSON-API Error Object](https://jsonapi.org/format/#error-objects) structure. The exception classes can be extended to pre-populate fields.

## Services

Currently, the filters are using the [HTTP Context](https://docs.nestjs.com/fundamentals/execution-context#execution-context). This should be abstracting using the [Execution Context](https://docs.nestjs.com/fundamentals/execution-context#current-application-context) to allow the same filters to work with services.

## Running unit tests

Run `ng test shared-exceptions` to exe cute the unit tests via [Jest](https://jestjs.io).
These filters catch errors that are thrown by the mongo database. In general, these should all be returning 500, even if they are casued by bad input from the user. This is because we want the API validation layer to catch the user errors and provide useful information to the user there. If the bad input reaches the database, then we have missed something, and its a 500 error for the purpose of logging/alerts. 
The one exception is the key conflict error, since a 409 fits better here, and checking this at the api layer has the same effect.# shared-utils-routes

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-utils-routes` to execute the unit tests.
# shared-services

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-services` to execute the unit tests.
# shared-angular

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-angular` to execute the unit tests.
# shared-content

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-content` to execute the unit tests.
# shared-ui

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-ui` to execute the unit tests.
# shared-error-handler

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-error-handler` to execute the unit tests.
# shared-nats-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test shared-nats-client` to execute the unit tests via [Jest](https://jestjs.io).
# shared-pwa

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-pwa` to execute the unit tests.
# shared-assets
# shared-debug

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-debug` to execute the unit tests.
# shared-environments

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test shared-environments` to execute the unit tests.
# shared-styles
# account-management

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test account-management` to execute the unit tests via [Jest](https://jestjs.io).
# api-angular-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test api-angular-client` to execute the unit tests.
# api-nest-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test api-nest-client` to execute the unit tests via [Jest](https://jestjs.io).
# hsds-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test hsds-client` to execute the unit tests via [Jest](https://jestjs.io).
# simulation-runs-service

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test simulation-runs-service` to execute the unit tests.
# simulation-runs-ui

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test simulation-runs-ui` to execute the unit tests.
# simulation-runs-viz

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test simulation-runs-viz` to execute the unit tests.
# hdf5-api-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test hdf5-api-client` to execute the unit tests via [Jest](https://jestjs.io).
## @

### Building

To install the required dependencies and to build the typescript sources run:
```
npm install
npm run build
```

#### General usage

In your Nestjs project:


```
// without configuring providers
import { ApiModule } from '';
import { HttpModule } from '@nestjs/common';

@Module({
    imports: [
        ApiModule,
        HttpModule
    ],
    providers: []
})
export class AppModule {}
```

```
// configuring providers
import { ApiModule, Configuration, ConfigurationParameters } from '';

export function apiConfigFactory (): Configuration => {
  const params: ConfigurationParameters = {
    // set configuration parameters here.
  }
  return new Configuration(params);
}

@Module({
    imports: [ ApiModule.forRoot(apiConfigFactory) ],
    declarations: [ AppComponent ],
    providers: [],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```

```
import { DefaultApi } from '';

export class AppComponent {
    constructor(private apiGateway: DefaultApi) { }
}
```

Note: The ApiModule a dynamic module and instantiated once app wide.
This is to ensure that all services are treated as singletons.

#### Using multiple swagger files / APIs / ApiModules
In order to use multiple `ApiModules` generated from different swagger files,
you can create an alias name when importing the modules
in order to avoid naming conflicts:
```
import { ApiModule } from 'my-api-path';
import { ApiModule as OtherApiModule } from 'my-other-api-path';
import { HttpModule } from '@nestjs/common';

@Module({
  imports: [
    ApiModule,
    OtherApiModule,
    HttpModule
  ]
})
export class AppModule {

}
```


### Set service base path
If different than the generated base path, during app bootstrap, you can provide the base path to your service. 

```
import { BASE_PATH } from '';

bootstrap(AppComponent, [
    { provide: BASE_PATH, useValue: 'https://your-web-service.com' },
]);
```
or

```
import { BASE_PATH } from '';

@Module({
    imports: [],
    declarations: [ AppComponent ],
    providers: [ provide: BASE_PATH, useValue: 'https://your-web-service.com' ],
    bootstrap: [ AppComponent ]
})
export class AppModule {}
```


#### Using @nestjs/cli
First extend your `src/environments/*.ts` files by adding the corresponding base path:

```
export const environment = {
  production: false,
  API_BASE_PATH: 'http://127.0.0.1:8080'
};
```

In the src/app/app.module.ts:
```
import { BASE_PATH } from '';
import { environment } from '../environments/environment';

@Module({
  declarations: [
    AppComponent
  ],
  imports: [ ],
  providers: [
    { 
      provide: 'BASE_PATH', 
      useValue: environment.API_BASE_PATH 
    }
  ]
})
export class AppModule { }
```
# messages

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test messages` to execute the unit tests via [Jest](https://jestjs.io).
# mail-service-client

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `ng test mail-service-client` to execute the unit tests via [Jest](https://jestjs.io).
# analytics-angular-analytics

This library was generated with [Nx](https://nx.dev).

## Running unit tests

Run `nx test analytics-angular-analytics` to execute the unit tests.
# S3 Bucket Configuration

## Development Bucket

## Production Bucket

## Temp Development Bucket

## Temp Prod Bucket

## CORS settings
The cors settings are used to configure the CORS policy for the bucket. The cors policy is contained in the `cors.py` file.

## Lifecycle Rules

The configuration for the retention policies and the lifecycle management is contained in the retention.py file.

### Generated COMBINE Archives

Generated combine archives are created to help users get valid omex archvies for use in runBiosimulations. The users upload a modelfile and then provide some information through a webform. We will delete these objects after 1 day to prevent any accidental retention of personal information.
If the omex archive is then used to run a simulation/create some resource, a copy will be saved through that process according to those policies.
**What new features does this PR implement?**
Please summarize the features that this PR implements. As relevant, please indicate the issues that the PR closes.

* Adds ABC (closes #X)
* Adds XYZ (closes #X)

**What bugs does this PR fix?**
Please summarize the bugs that this PR fixes. As relevant, please indicate the issues that the PR closes.

* Fixes ABC (closes #X)
* Fixes XYZ (closes #X)

**Does this PR introduce any additional changes?**
Please summarize any additional changes that this PR introduces.

* Corrects typos in documentation
* Changes theme of documentation

**How have you tested this PR?**
Please summarize the tests that you implemented to test this PR.

* Added test for ABC
* Added test for XYZ

**Additional information**
Please describe any information needed to review this PR.
---
title: Example simulation runs need to be fixed
labels: bug, dispatch
---

{{ env.OUTPUT }}
<!-- See https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples%3A-advice --> <!-- markdownlint-disable MD033 MD041 -->
<details><summary>If you see a bunch of garbage</summary>

If it relates to a ...
<details><summary>well-formed pattern</summary>

See if there's a [pattern](https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples:-patterns) that would match it.

If not, try writing one and adding it to the `patterns.txt` file.

Patterns are Perl 5 Regular Expressions - you can [test](
https://www.regexplanet.com/advanced/perl/) yours before committing to verify it will match your lines.

Note that patterns can't match multiline strings.
</details>
<details><summary>binary-ish string</summary>

Please add a file path to the `excludes.txt` file instead of just accepting the garbage.

File paths are Perl 5 Regular Expressions - you can [test](
https://www.regexplanet.com/advanced/perl/) yours before committing to verify it will match your files.

`^` refers to the file's path from the root of the repository, so `^README\.md$` would exclude [README.md](
../tree/HEAD/README.md) (on whichever branch you're using).
</details>

</details>
# check-spelling/check-spelling configuration

File | Purpose | Format | Info
-|-|-|-
[dictionary.txt](dictionary.txt) | Replacement dictionary (creating this file will override the default dictionary) | one word per line | [dictionary](https://github.com/check-spelling/check-spelling/wiki/Configuration#dictionary)
[allow.txt](allow.txt) | Add words to the dictionary | one word per line (only letters and `'`s allowed) | [allow](https://github.com/check-spelling/check-spelling/wiki/Configuration#allow)
[reject.txt](reject.txt) | Remove words from the dictionary (after allow) | grep pattern matching whole dictionary words | [reject](https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples%3A-reject)
[excludes.txt](excludes.txt) | Files to ignore entirely | perl regular expression | [excludes](https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples%3A-excludes)
[only.txt](only.txt) | Only check matching files (applied after excludes) | perl regular expression | [only](https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples%3A-only)
[patterns.txt](patterns.txt) | Patterns to ignore from checked lines | perl regular expression (order matters, first match wins) | [patterns](https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples%3A-patterns)
[expect.txt](expect.txt) | Expected words that aren't in the dictionary | one word per line (sorted, alphabetically) | [expect](https://github.com/check-spelling/check-spelling/wiki/Configuration#expect)
[advice.md](advice.md) | Supplement for GitHub comment when unrecognized words are found | GitHub Markdown | [advice](https://github.com/check-spelling/check-spelling/wiki/Configuration-Examples%3A-advice)

Note: you can replace any of these files with a directory by the same name (minus the suffix)
and then include multiple files inside that directory (with that suffix) to merge multiple files together.
# COMBINE-API 

COMBINE-API is an HTTP API for working with COMBINE/OMEX archives and other COMBINE formats, such as the OMEX manifest, OMEX metadata, and SED-ML formats and several model formats such as BNGL, CellML, LEMS, NeuroML, RBA XML, SBML, Smoldyn, and XPP. This API is specified using OpenAPI and implemented in Python using [connexion](https://github.com/zalando/connexion). This API is currently deployed at https://combine.api.biosimulations.org and https://combine.api.biosimulations.dev.

## Editing the specifications of the API
The [specifications](src/spec/spec.yml) of the API were developed using [Apicurio studio](https://www.apicur.io/studio/), a free web-based application for building [OpenAPI](https://swagger.io/specification/)-compliant specifications for HTTP APIs. The specifications of the API should be edited using Apicurio studio, exported from Apicurio in YAML format, and saved to [`src/spec/spec.yml`](src/spec/spec.yml) in this repository. The Apicurio project for the API is owned by Bilal Shaikh. Contact Bilal for priviledges to edit the project.

Please follow these steps to edit the specifications of the API:
1. Login to [Apicurio studio](https://www.apicur.io/studio/).
2. Use Apicurio studio to edit the API as needed.
3. Navigate to the [Apicurio landing page page for the API](https://studio.apicur.io/apis/45980).
4. Click the '...' button at the top-right of the box for the API.
5. Click the 'Download' menu option.
6. Select 'Format' 'YAML'.
7. Select 'References' 'Deference All External $refs'.
8. Click the 'Download' button.
9. Save the exported YAML file to [`src/spec/spec.yml`](src/spec/spec.yml) in this repository.
10. Commit the modified specifications file to this repository.

## Editing the implementation of the API
The API is implemented in Python using [connexion](https://github.com/zalando/connexion). To edit the operation for a path, open the Python function for the path in your favorite Python editor and edit the function. The specifications of the API indicate the Python functions which handle each path (`operationId` of each HTTP method of each path). For example, the `operationId` for the `/health` path is `src.handlers.health.handler`, which corresponds to the `handler` function in the [`src.handlers.health`](src/handlers/health.py) Python module. connexion automatically maps each path to the functions indicated by the `operationId` attributes in the specifications of the API.

## Installation and execution
As described below, the COMBINE API is deployed as a Docker image. The [`Dockerfile`](Dockerfile) for this image is the authoritative description of how to install and execute the API. The `Dockerfile` uses [pipenv](https://pipenv.pypa.io/) to install the required Python packages outlined in the [`Pipfile`](Dockerfile-assets/Pipfile) and [`Pipfile.lock`](Dockerfile-assets/Pipfile.lock) files. Note, this `Pipfile` does not completely describe the requirements for the API. The `Dockerfile` includes additional operations which cannot be achieved using pipenv because some of the required Python packages for the API require additional OS packages and because bugs in pipenv currently prevent pipenv from installing one of the required Python packages for the API.

### Recommended local installation
Below is an outline of how to install the API into a local Python environment managed with pipenv in a Linux machine.

1. Change to this directory (`apps/combine-api` subdirectory of this repository).
2. Install the OS packages outlined in `Dockerfile`:
    ```bash
    apt-get install default-jre perl r-base ...
    ```
3. Use pipenv to create a Python environment with the required Python packages:
    ```bash
    apt-get install python3 python3-pip
    pip3 install pipenv
    cd Dockerfile-assets
    pipenv install
    ```
4. Install the additional Python packages outlined at the end of the `Dockerfile` (steps after the `pipenv install ...` step). These packages must be installed separately from pipenv because a bug in pipenv currently prevents pipenv from installing them.

### Adding additional requirements for the API
1. Follow the recommended installation steps above.
2. Change to this directory (`apps/combine-api` subdirectory of this repository).
3. Use `apt` or other methods to install any additional required OS packages.
4. Add these additional required OS packages to `Dockerfile`.
5. Change to the directory for the Python environment for the API (`apps/combine-api/Dockerfile-assets` subdirectory of this repository).
6. Run `pipenv install ...` to install any additional required Python packages.
7. Commit the modified `Dockerfile`, `Pipefile`, and `Pipfile.lock` files to this repository.

### Recommended execution of a local development server
1. Follow the recommended installation steps above.
2. Change to the directory for the Python environment for the API (`apps/combine-api/Dockerfile-assets` subdirectory of this repository).
3. Activate this Python environment by running `pipenv shell`.
4. Change to the root directory for this API by running (`cd ..`).
5. Execute a server for the API by running `python -m src`.
6. Navigate your browser to the URL printed to your console (the default is currently http://0.0.0.0:3333/).

### Building a Docker image for the API
1. Change to this directory (`apps/combine-api` subdirectory of this repository).
2. Run `docker build --tag ghcr.io/biosimulations/combine-api .`.

### Pulling a published Docker image for the API
As described below, Docker images for the API are available from the GitHub Container Registry ([`ghcr.io/biosimulations/combine-api`](https://github.com/biosimulations/biosimulations/pkgs/container/combine-api)). These images can be pulled by running `docker pull ghcr.io/biosimulations/combine-api`.

### Executing a Docker image for the API
Follow these steps to run an image for the API:
1. Execute the image on local port `3333` by running `docker run -it --rm -p 127.0.0.1:3333:3333 ghcr.io/biosimulations/combine-api`. The second `3333` must match the port which the API is running at in the image (the default is currently `3333`).
2. Navigate your browser to `http://127.0.0.1:3333`.

## Linting the API
The API can be linted using [flake8](https://flake8.pycqa.org/en/latest/). 

### Linting the API locally
1. Follow the recommended installation steps above.
2. Change to the directory for the Python environment for the API (`apps/combine-api/Dockerfile-assets` subdirectory of this repository).
3. Activate this Python environment by running `pipenv shell`.
4. Change to the root directory for this repository (`cd ../../..`).
5. Run `flake8`.

### Linting execution by the CI system
The CI system lints the API by running `nx run combine-api:lint` from the root directory of this repository.

## Testing the API
The tests for the API are located in the [`tests`](tests) subdirectory of this directory. The tests can be executed using [pytest](https://docs.pytest.org/). Coverage can be accessed using [pytest-cov](https://pytest-cov.readthedocs.io/). 

### Testing the API locally
1. Follow the recommended installation steps above.
2. Change to the directory for the Python environment for the API (`apps/combine-api/Dockerfile-assets` subdirectory of this repository).
3. Activate this Python environment by running `pipenv shell`.
4. Change to the root directory for this API by running (`cd ..`).
5. Run `python -m pytest tests/`.
   a. To measure the coverage of the tests, add the `--cov src` option.
   b. To compile the coverage report to HTML, run `coverage html`.
   c. To view the coverage report, navigate your browser to the `htmlcov/` subdirectory of the directory for the Python environment for the API (`apps/combine-api/Dockerfile-assets/htmlcov/index.html`).

### Test execution by the CI system
The CI system executes these tests by running `nx run combine-api:test` from the root directory of this repository.

## Deployment
The COMBINE API is deployed as a Docker image. [`Dockerfile`](Dockerfile) is the Dockerfile for this image. This Dockerfile uses several files in the [`Dockerfile-assets`](Dockerfile-assets) subdirectory of this directory. This includes [pipenv](https://pipenv.pypa.io/) [`Pipfile`](Dockerfile-assets/Pipfile) and [`Pipfile.lock`](Dockerfile-assets/Pipfile.lock) files which describe most of the Python packages required for the API. Note, this `Pipfile` does not completely describe the requirements for the API. The `Dockerfile` includes additional operations which cannot be achieved using pipenv because some of the required Python packages for the API require additional OS packages and because bugs in pipenv currently prevent pipenv from installing one of the required Python packages for the API.

The CI/CD system builds this image by running `nx docker combine-api` from the root directory of this repository and pushes it to the GitHub Container Registry.
# BioSimulations and BioSimulators documentation

## Motivation and goals

More comprehensive and more predictive models have the potential to advance biology, bioengineering, and medicine. Building more predictive models will likely require the collaborative efforts of many investigators. This requires teams to be able to share and reuse models, simulations, and simulation tools. Despite extensive efforts to develop standards formats such as the [Systems Biology Markup Language (SBML)](http://sbml.org/), the [Simulation Experiment Description Markup Language (SED-ML)](https://www.sed-ml.org/), and the [COMBINE/OMEX archive format](https://combinearchive.org/); ontologies such as the [Kinetic Simulation Algorithm Ontology (KiSAO)](https://github.com/SED-ML/KiSAO/) and the [Systems Biology Ontology (SBO)](https://github.com/EBI-BioModels/SBO); and repositories such as [BioModels](http://biomodels.net/) and the [Physiome Model Repository](https://models.physiomeproject.org/), it is still often difficult to share and reuse models and simulations. One challenge to sharing and reusing models is the disparate modeling formalisms, simulation algorithms, model formats, simulation tools, and repositories for different types of models and biological systems. This proliferation of methods, formats, repositories, and tools makes it difficult, especially for non-experts, to find models and simulation tools for reusing them. In addition, the existing model repositories have limited capabilities for sharing associated resources such as training data, simulation experiments, and interactive visualizations for their results.

## Key features

BioSimulations and BioSimulators address these challenges by making it easier for researchers to share and reuse simulations and simulation tools. 

* Central portals for publishing and discovering simulation projects and tools. BioSimulations and BioSimulators provide central portals for sharing and discovering models, simulations, data visualizations, and simulation tools across a broad range of modeling frameworks, model formats, simulation algorithms, simulation tools, and data visualizations.
* Guides to finding simulation tools: BioSimulators provides investigators centralized information about the capabilities of simulation tools. This information can be used to find appropriate tools for specific projects. In addition, runBioSimulations provides services for automatically recommending tools for specific COMBINE/OMEX archives, modeling frameworks, or simulation algorithms.
* Web-based tools for running simulations and interactively exploring their results. runBioSimulations provides a simple interface for executing simulations and interactively visualizing their results. runBioSimulations also supports a broad range of modeling methods and formats.
* Simple tools for reusing simulation projects. runBioSimulations provides a simple web interface for modifying simulation experiments and running new simulations.
* Transparent simulation. By building upon the COMBINE/OMEX archive format, KiSAO, and SED-ML, the details of each simulation experiment are fully transparent. This helps investigators understand and reproduce simulation experiments.
* Seamless integration of local development and publication. By building on containerization, runBioSimulations provides authors the ability to preview their simulations using their own computers prior to publication. Similarly, this enables investigators to use their own resources to explore simulations obtained from BioSimulations. In particular, this avoids duplicate effort in using simulations in multiple different environments and enables investigators to interactively debug problems using their own machines.

## Supported simulation methods

The BioSimulators platform supports a broad range of modeling frameworks, model formats, simulation algorithms, and simulation tools. Currently, the BioSimulators registry includes standardized interfaces to simulation tools for constraint-based (Flux Balance Analysis (FBA) and Resource Balance Analysis (RBA)), continuous kinetic (ordinary differential equations (ODE) and differential-algebraic equations (DAE)), discrete kinetic (e.g., Stochastic Simulation Algorithms (SSA)), logical, spatial, particle-based and hybrid models that are described in using several languages including the [BioNetGen Language (BNGL)](https://bionetgen.org), [CellML](https://cellml.org), the [GINsim](http://ginsim.org/) Markup Language, [NeuroML](https://neuroml.org/)/[Low Entropy Model Specification Langauge (LEMS)](https://lems.github.io/LEMS/), the [RBA XML format](https://sysbioinra.github.io/RBApy/), the [Systems Biology Markup Language (SBML)](https://sbml.org) including the Flux Balance Constraints and Qualitative Models Packages, the [Smoldyn](http://www.smoldyn.org/) simulation configuration format, and the XPP [ODE](http://www.math.pitt.edu/~bard/xpp/help/xppodes.html) format. This encompasses over 60 simulation algorithm algorithms with over 20 simulation tools. More information about the available simulation methods available through BioSimulators is available at [https://biosimulators.org](https://biosimulators.org). These simulation capabilities are available through runBioSimulations and BioSimulations. Further, the community can extend BioSimulators' capabilities by contributing additional simulation tools. More information, tutorials, and examples are available from BioSimulators and in this documentation.

<div class="logos">
<div class="logos-row">
    <a href="https://www.bionetgen.org" rel="noopener" target="_blank" title="BNGL">
    <img
        class="zoom"
        src="/assets/images/about/partners/bionetgen.png"
    />
    </a>

    <a href="https://www.cellml.org/" rel="noopener" target="_blank" title="CellML">
    <img class="zoom" src="/assets/images/about/partners/cellml.svg" />
    </a>

    <a href="http://ginsim.org/" rel="noopener" target="_blank" title="GINsim">
    <img class="zoom" src="/assets/images/about/partners/ginsim.svg" />
    </a>

    <a href="https://neuroml.org/" rel="noopener" target="_blank" title="NeuroML">
    <img class="zoom" src="/assets/images/about/partners/neuroml.svg" />
    </a>

    <!--
    <a href="http://www.pharmml.org/" rel="noopener" target="_blank" title="pharmML">
    <img class="zoom" src="/assets/images/about/partners/pharmml.svg" />
    </a>
    -->

    <a href="https://rba.inrae.fr" rel="noopener" target="_blank" title="RBA">
    <img class="zoom" src="/assets/images/about/partners/rba.png" />
    </a>

    <a href="http://sbml.org" rel="noopener" target="_blank" title="SBML">
    <img class="zoom" src="/assets/images/about/partners/sbml.svg" />
    </a>
</div>

<div class="logos-row">
    <a href="http://www.ebi.ac.uk/sbo/" rel="noopener" target="_blank" title="SBO">
    <img class="zoom" src="/assets/images/about/partners/sbo.png" />
    </a>

    <a href="https://sed-ml.org/" rel="noopener" target="_blank" title="SED-ML">
    <img class="zoom" src="/assets/images/about/partners/sed-ml.svg" />
    </a>

    <a
    href="http://co.mbine.org/standards/kisao"
    rel="noopener" target="_blank"
    title="KiSAO"
    >
    <img class="zoom" src="/assets/images/about/partners/kisao.svg" />
    </a>

    <a href="https://escher.github.io/" rel="noopener" target="_blank" title="Escher">
    <img class="zoom" src="/assets/images/about/partners/escher.svg" />
    </a>

    <a href="https://sbgn.github.io/" rel="noopener" target="_blank" title="SBGN">
    <img class="zoom" src="/assets/images/about/partners/sbgn.png" />
    </a>

    <a href="https://vega.github.io/vega/" rel="noopener" target="_blank" title="Vega">
    <img class="zoom" src="/assets/images/about/partners/vega.svg" />
    </a>

    <a
    href="https://co.mbine.org/standards/omex"
    rel="noopener" target="_blank"
    title="OMEX"
    >
    <img class="zoom" src="/assets/images/about/partners/omex.svg" />
    </a>
</div>

<div class="logos-row">
    <a
    href="http://amici-dev.github.io/AMICI/"
    rel="noopener" target="_blank"
    title="AMICI"
    >
    <img class="zoom" src="/assets/images/about/partners/amici.svg" />
    </a>

    <a href="https://bionetgen.org" rel="noopener" target="_blank" title="BioNetGen">
    <img
        class="zoom"
        src="/assets/images/about/partners/bionetgen.png"
    />
    </a>

    <!--
    <a
    href="https://cayenne.readthedocs.io/"
    rel="noopener" target="_blank"
    title="Cayenne"
    >
    <img class="zoom" src="/assets/images/about/partners/cayenne.png" />
    </a>
    -->

    <a
    href="https://opencobra.github.io/cobrapy/"
    rel="noopener" target="_blank"
    title="COBRApy"
    >
    <img class="zoom" src="/assets/images/about/partners/cobrapy.svg" />
    </a>

    <a href="http://copasi.org/" rel="noopener" target="_blank" title="COPASI">
    <img class="zoom" src="/assets/images/about/partners/copasi.svg" />
    </a>

    <a
    href="https://gillespy2.github.io/GillesPy2/"
    rel="noopener" target="_blank"
    title="GillesPy2"
    >
    <img
        class="zoom"
        src="/assets/images/about/partners/gillespy2.svg"
    />
    </a>

    <!--
    <a
    href="https://github.com/MyersResearchGroup/iBioSim"
    rel="noopener" target="_blank"
    title="iBioSim"
    >
    <img class="zoom" src="/assets/images/about/partners/ibiosim.svg" />
    </a>
    -->

    <a
    href="https://masspy.readthedocs.io/"
    rel="noopener" target="_blank"
    title="MASSpy"
    >
    <img class="zoom" src="/assets/images/about/partners/masspy.svg" />
    </a>

    <a href="http://www.netpyne.org/" rel="noopener" target="_blank" title="NetPyNe">
    <img class="zoom" src="/assets/images/about/partners/netpyne.png" />
    </a>
</div>

<div class="logos-row">
    <a
    href="https://sysbioinra.github.io/RBApy/"
    rel="noopener" target="_blank"
    title="RBApy"
    >
    <img class="zoom" src="/assets/images/about/partners/rbapy.svg" />
    </a>

    <a
    href="http://pysces.sourceforge.net/"
    rel="noopener" target="_blank"
    title="PySCeS"
    >
    <img class="zoom" src="/assets/images/about/partners/pysces.svg" />
    </a>

    <a
    href="http://tellurium.analogmachine.org/"
    rel="noopener" target="_blank"
    title="tellurium/libRoadRunner"
    >
    <img
        class="zoom"
        src="/assets/images/about/partners/libroadrunner.svg"
    />
    </a>

    <a href="https://vcell.org/" rel="noopener" target="_blank" title="VCell">
    <img class="zoom" src="/assets/images/about/partners/vcell.svg" />
</a>
</div>
</div>
# Changelog

## [8.8.0](https://github.com/biosimulations/biosimulations/compare/v8.7.1...v8.8.0) (2022-01-08)


### Features

* added workflow to delete temporary COMBINE archives ([0d722e9](https://github.com/biosimulations/biosimulations/commit/0d722e939b19726c0f91c310d8a0a6a31217aed8))
* **api,combine-api,dispatch,platform:** added support for references for projects ([a544969](https://github.com/biosimulations/biosimulations/commit/a5449691805801726b643285a018dbff77e06a00))

## [8.7.1](https://github.com/biosimulations/biosimulations/compare/v8.7.0...v8.7.1) (2022-01-06)


### Bug Fixes

* add script ignore false to sharp install ([9b90d9b](https://github.com/biosimulations/biosimulations/commit/9b90d9bfeefe1e8b1044e9e2685e5de15bb3b2e4))
* **dispatch-service:** added dependency for sharp to Dockerfile ([ba344c7](https://github.com/biosimulations/biosimulations/commit/ba344c7649c97b842b9794e3f58bd8af3cbe2712))
* **dispatch,platform:** fixed name and URL for log format in files tab ([f13647f](https://github.com/biosimulations/biosimulations/commit/f13647fcaac05b5abdfb5ad4a754c3f4c18586ac))

## [8.7.0](https://github.com/biosimulations/biosimulations/compare/v8.6.0...v8.7.0) (2022-01-06)


### Bug Fixes

* **dispatch-service:** separated input and output files for simulation runs ([dc98a46](https://github.com/biosimulations/biosimulations/commit/dc98a4604add2541e1214afa4701eb829c634a3f))
* **dispatch,platform,ui:** corrected layout of metadata columns ([8044bc6](https://github.com/biosimulations/biosimulations/commit/8044bc62d06057696b1a12d22d61a9287d1dbd7f))
* **dispatch,platform:** added modeling methods to metadata ([6da720e](https://github.com/biosimulations/biosimulations/commit/6da720ef4504515f9d9eafdb11e8a76d56a448cc))
* **dispatch:** add rel noopenor for external links ([9efa059](https://github.com/biosimulations/biosimulations/commit/9efa059a50435583859f17ff05476694a84a5e59))
* **ui:** fixed setting of open control panel in table controls ([ca529a3](https://github.com/biosimulations/biosimulations/commit/ca529a3be58ca72ab4791732de69ff0efa8383fa))


### Features

* **combine-api:** added utility methods for reading S3 files ([4e3df4f](https://github.com/biosimulations/biosimulations/commit/4e3df4fee4247f1d4ae371f634cb6af8084b93ee))
* **dispatch-service:** added support for SLURM constraints ([df2acba](https://github.com/biosimulations/biosimulations/commit/df2acba2839c4def3139dc0a889af4242889fbde))
* **dispatch,platform,ui:** added markdown rendering for project descriptions ([818d267](https://github.com/biosimulations/biosimulations/commit/818d2673f20c1a87f235618cae242c129b49630f))
* **dispatch,platform,ui:** interleaved metadata about files into files tab ([9807406](https://github.com/biosimulations/biosimulations/commit/9807406cbcfc1c78b793d2068a98b4693725741e))
* **dispatch:** made it easier to get errors with simulation projects ([4c7635f](https://github.com/biosimulations/biosimulations/commit/4c7635fb9997d8cb33f4c2051bb6ffeb9958c42e))
* **ontology:** added additional formats used by Physiome ([3e951d0](https://github.com/biosimulations/biosimulations/commit/3e951d02cd6633aca35d444ce60dd076dbcb1dc8))


### Performance Improvements

* **api:** added caching for getting properties of ontology terms ([c273e88](https://github.com/biosimulations/biosimulations/commit/c273e88cb31efeaf7872aa48ad1aa38a05524b80))
* **dispatch,platform:** reduced thumbnail image sizes ([049af36](https://github.com/biosimulations/biosimulations/commit/049af36affaf05b56a35834c0c342803d74ba46c))

## [8.6.0](https://github.com/biosimulations/biosimulations/compare/v8.5.6...v8.6.0) (2021-12-31)


### Features

* **ontology:** added formats used by Physiome model repository ([b407218](https://github.com/biosimulations/biosimulations/commit/b40721893c8592e1418f21fe78462d1978270056))

## [8.5.6](https://github.com/biosimulations/biosimulations/compare/v8.5.5...v8.5.6) (2021-12-27)


### Bug Fixes

* update package lock version ([b91c669](https://github.com/biosimulations/biosimulations/commit/b91c669745abee0599cfd002fd18f557acd8e155))

## [8.5.5](https://github.com/biosimulations/biosimulations/compare/v8.5.4...v8.5.5) (2021-12-23)


### Bug Fixes

* **api:** tried to correct logging for errors in uploading archives ([b492f02](https://github.com/biosimulations/biosimulations/commit/b492f02c04da3bc6bf1bdd7ec022ce85187ddfa1))

## [8.5.4](https://github.com/biosimulations/biosimulations/compare/v8.5.3...v8.5.4) (2021-12-23)


### Bug Fixes

* improved error logging, increasing AWS timeout ([cccf802](https://github.com/biosimulations/biosimulations/commit/cccf802c76f154ed09dc53a7b1a41df286d935a8))

## [8.5.3](https://github.com/biosimulations/biosimulations/compare/v8.5.2...v8.5.3) (2021-12-23)


### Bug Fixes

* **dispatch-service:** fixed retrying of files and specs to avoid conflicts on incomplete posts ([620f0eb](https://github.com/biosimulations/biosimulations/commit/620f0eb59a1f347ea9bbbaaf0b3d0276654d8611))

## [8.5.2](https://github.com/biosimulations/biosimulations/compare/v8.5.1...v8.5.2) (2021-12-23)


### Bug Fixes

* **api,dispatch-service:** fixed logging for complete processor failures ([f2119d4](https://github.com/biosimulations/biosimulations/commit/f2119d46741eb85f89a91d84eeee93903900acc8))

## [8.5.1](https://github.com/biosimulations/biosimulations/compare/v8.5.0...v8.5.1) (2021-12-23)


### Bug Fixes

* **dispatch-service:** retry project publication after run completion ([e686664](https://github.com/biosimulations/biosimulations/commit/e68666451fc61bd27dd25e3775a4477fca5385e1))

## [8.5.0](https://github.com/biosimulations/biosimulations/compare/v8.4.1...v8.5.0) (2021-12-22)


### Bug Fixes

* **simulators-api:** fixed updating of updated timestamp; addresses [#3878](https://github.com/biosimulations/biosimulations/issues/3878) ([13ee425](https://github.com/biosimulations/biosimulations/commit/13ee42542af3f15c13954e7873313581cd76647c))
* **ui,dispatch:** fixed unselecting files; closes [#3875](https://github.com/biosimulations/biosimulations/issues/3875) ([2f84e18](https://github.com/biosimulations/biosimulations/commit/2f84e187ee1ea2c865925d486b5235682814e0da))


### Features

* **ui:** added instructions to refresh on failures ([0181c3b](https://github.com/biosimulations/biosimulations/commit/0181c3be83cfe5602d258ca8276599348c123f9e))

## [8.4.1](https://github.com/biosimulations/biosimulations/compare/v8.4.0...v8.4.1) (2021-12-22)


### Bug Fixes

* **combine-api:** fixed reading shared configuration ([5a83e6e](https://github.com/biosimulations/biosimulations/commit/5a83e6e0f4802eed023781dac8881f4c8e251a2e))

## [8.4.0](https://github.com/biosimulations/biosimulations/compare/v8.3.0...v8.4.0) (2021-12-22)


### Bug Fixes

* added project id, owner to CompleteJob for failures ([0a8d834](https://github.com/biosimulations/biosimulations/commit/0a8d834e100fc396da08403e5cef6e8c703c133c))
* **api,dispatch:** fixed data model for simulation results ([7cd2b5e](https://github.com/biosimulations/biosimulations/commit/7cd2b5e2a3d37e7544c16b383ee6594b65afdf69))
* **api,dispatch:** fixed data type for simulation results ([610cbc9](https://github.com/biosimulations/biosimulations/commit/610cbc97530b8b9965a98a356beb23e27fbf3cdc))
* **api:** changed file URL validation to allow un-encoded URLs ([f78bc43](https://github.com/biosimulations/biosimulations/commit/f78bc43fd89bb990841786eacb5118b5645628c4))
* fixed external simulators API endpoint ([25b3a3c](https://github.com/biosimulations/biosimulations/commit/25b3a3c0aabbf14949cf5386438b4bee87c95a03))
* fixed spinner for loading table data ([e2e1314](https://github.com/biosimulations/biosimulations/commit/e2e13149988f7cd74f95849f9c473c914d44f5ee))


### Features

* added checks that S3 files were deleted ([390300c](https://github.com/biosimulations/biosimulations/commit/390300c6039cede0176f7a5824b01687dbae5a94))
* **api:** added cache for project summaries ([6fc5bb7](https://github.com/biosimulations/biosimulations/commit/6fc5bb705a3cadb3ee2379145b3a362f8abac354))
* **api:** added checks that S3 files were deleted ([5d3ccb9](https://github.com/biosimulations/biosimulations/commit/5d3ccb91f964b11e9ffca1db1cb5f51d2c1d4389))
* **api:** added handled for NaN and Inf from HSDS ([7802310](https://github.com/biosimulations/biosimulations/commit/78023109c8a0b3c1c6d4e80ec65aec680b8a6172))
* **combine-api:** relaxed required metadata for simulation projects ([2be28f1](https://github.com/biosimulations/biosimulations/commit/2be28f1a896256e29ea7f9387e66f6283538e094))
* directed href targets ([d609f80](https://github.com/biosimulations/biosimulations/commit/d609f80c71549213eaac6b9b1fd8a43b7d2c038c))

## [8.3.0](https://github.com/biosimulations/biosimulations/compare/v8.2.1...v8.3.0) (2021-12-20)


### Bug Fixes

* **dispatch-service:** fix dispatch service post limit ([95653f8](https://github.com/biosimulations/biosimulations/commit/95653f88710f37babc19ccfa85b7c7caabc9fec3)), closes [#3828](https://github.com/biosimulations/biosimulations/issues/3828)
* **dispatch:** correct security issue with untrusted html input ([8b464b6](https://github.com/biosimulations/biosimulations/commit/8b464b64ecc1992031e49b76a4765ce68a70b7f4))
* **ui:** fixed spinner exit for table component ([f26802c](https://github.com/biosimulations/biosimulations/commit/f26802cf284df664509f2052a5f29a4443ee9a36))


### Features

* **api:** added project summary caching at creation and updating ([20d99f2](https://github.com/biosimulations/biosimulations/commit/20d99f2bcccba153a4aae10b9872c212935563f6))
* **dispatch:** added example simulation runs for Brian 2 ([710df68](https://github.com/biosimulations/biosimulations/commit/710df689dc79f3ac837ddda2034398371eb5f084))

## [8.2.1](https://github.com/biosimulations/biosimulations/compare/v8.2.0...v8.2.1) (2021-12-16)


### Bug Fixes

* add gtag snippet to dispatch and simulators ([f1b6332](https://github.com/biosimulations/biosimulations/commit/f1b633298b064a3a9d6ce8bdc404fd815ead5a5c))
* **config:** add default server limit to config ([496b430](https://github.com/biosimulations/biosimulations/commit/496b430efd96e7d2b13b102ea8f7ef9d25b8e35a)), closes [#3828](https://github.com/biosimulations/biosimulations/issues/3828)

## [8.2.0](https://github.com/biosimulations/biosimulations/compare/v8.1.0...v8.2.0) (2021-12-16)


### Bug Fixes

* fixed log validation ([86fc30d](https://github.com/biosimulations/biosimulations/commit/86fc30d388cf1ee170952be7efabdbd5bc5faca7))


### Features

* add angular analytics package ([3363f7b](https://github.com/biosimulations/biosimulations/commit/3363f7bcae10fcc07383e36c617bb960b054f380))
* added implementation of analytics and user consent ([2d87bb1](https://github.com/biosimulations/biosimulations/commit/2d87bb16abe21c9af02f402e3cfdb9265b3605e6))
* **dispatch,platform,simulators:** add cookie consent and privacy settings to frontend apps ([e84cdea](https://github.com/biosimulations/biosimulations/commit/e84cdeaa8a230a068afbb490dc33a796a441cc59))

## [8.1.0](https://github.com/biosimulations/biosimulations/compare/v8.0.0...v8.1.0) (2021-12-16)


### Bug Fixes

* fixed display of files in subdirectories ([ee62bbe](https://github.com/biosimulations/biosimulations/commit/ee62bbebf4cdeafcc3ec24f20a16bbc1178a0a82))
* fixed file size extraction for empty files ([b3dce39](https://github.com/biosimulations/biosimulations/commit/b3dce39fc83d9d846d1e1a87e5435e9cdc5078b2))
* fixed simulators view endpoint method ([e40f4cd](https://github.com/biosimulations/biosimulations/commit/e40f4cd29497f096d7341fd3c7b9c589b27f0093))
* set minimum time step to 1 ([e92ffa5](https://github.com/biosimulations/biosimulations/commit/e92ffa5de8a163b7e844a8f47cdd55e02b2b3155))


### Features

* **combine-api:** added error messages for invalid S3 bucket configuration ([626a2f2](https://github.com/biosimulations/biosimulations/commit/626a2f2a6deadaaa71517a09c19451f975316ad0))

## [8.0.0](https://github.com/biosimulations/biosimulations/compare/v7.0.0...v8.0.0) (2021-12-15)


### Bug Fixes

* **config:** change name of env variable to avoid clash ([caf69f4](https://github.com/biosimulations/biosimulations/commit/caf69f4ef11439ad6b6269ee41e904283219225c))
* **config:** correct endpoints for s3 contents path. Add some testing ([73f6034](https://github.com/biosimulations/biosimulations/commit/73f6034561fd88690a89f45d76dde9ed681b55ee)), closes [#3755](https://github.com/biosimulations/biosimulations/issues/3755)
* **config:** fix endpoints for ontology url ([f7ba9c5](https://github.com/biosimulations/biosimulations/commit/f7ba9c5f55048092280a0eeef0740bccc84d4c4e)), closes [#3771](https://github.com/biosimulations/biosimulations/issues/3771)
* **config:** fix prod file url ([5211afd](https://github.com/biosimulations/biosimulations/commit/5211afd21258dedd1483b7509270accfdc6f8dc8))
* **config:** fix storage health endpoint ([14d303c](https://github.com/biosimulations/biosimulations/commit/14d303c8c304468238b446acb1cbb7a4023e29b3))
* **dispatch-service:** corrected storage endpoint in sbatch sevice to external ([8eb34d4](https://github.com/biosimulations/biosimulations/commit/8eb34d4091e010b8c7b4ada5d29ceb1bf5c5a4d5))
* fixed broken links in documentation ([fe96ded](https://github.com/biosimulations/biosimulations/commit/fe96dedddd943974a5768a0f4439095f6bea8958))
* fixed broken links in documentation ([8de8a6f](https://github.com/biosimulations/biosimulations/commit/8de8a6f40ae046a95d5557d7d453220037eb9e27))
* removed invalid and unecessary workflow_call secret inputs ([b925dfe](https://github.com/biosimulations/biosimulations/commit/b925dfed75d24546ac5a01c6ac43fbf378b594b7))
* **simulators-api:** fix permissions for deletion of all simulation runs ([24aef29](https://github.com/biosimulations/biosimulations/commit/24aef2972670e4046707fb98072b0834f8435472)), closes [#3767](https://github.com/biosimulations/biosimulations/issues/3767)


### Code Refactoring

* **api:** simplify management of files ([e16c35f](https://github.com/biosimulations/biosimulations/commit/e16c35f374f5cbf0acda5f26288a1ed5c1ce04ec))


### Features

* **api,dispatch-service,dispatch:** expanded to full SED-ML data model ([2550def](https://github.com/biosimulations/biosimulations/commit/2550defc9918ae44dd2d5df53b6e2035ebcb7a00))
* **api:** improved error logging ([4d75193](https://github.com/biosimulations/biosimulations/commit/4d75193c68578f4dd48740795dfe9a630109366e))
* **auth:** add scope for deleting all simulation runs ([77a859c](https://github.com/biosimulations/biosimulations/commit/77a859caec20639adf77e6b0ac6d5ad00ecbb1ca))
* **combine-api:** expanded to full Python SED-ML data model ([912f8d5](https://github.com/biosimulations/biosimulations/commit/912f8d57bf205c76f7a12b4147dbf63e7bed0e89))
* **combine-api:** updated dependencies ([a7bf6e9](https://github.com/biosimulations/biosimulations/commit/a7bf6e9e4060658f62837e74172169abeb217f80))
* **config:** load endpoints dynamically if not in browser ([90ddc18](https://github.com/biosimulations/biosimulations/commit/90ddc1895fb71d04663c86965bf8fa3a6d1ec785)), closes [#3585](https://github.com/biosimulations/biosimulations/issues/3585)
* **dispatch-service:** added env variables to enable simulators to get number of CPUs ([fbdf8fa](https://github.com/biosimulations/biosimulations/commit/fbdf8fa354a6b06339d9831385c649a0e4d02e24))
* **dispatch,dispatch-service:** added support for multidimensional plots ([6681d25](https://github.com/biosimulations/biosimulations/commit/6681d25aed91b2bc51f5442ab929390765c9ed63))
* **dispatch,platform,api:** added filter and display of SED-ML file provenance ([d2a6c29](https://github.com/biosimulations/biosimulations/commit/d2a6c29a77ea1fc10289e576f10a6ef3bbb467b9))
* **dispatch:** extend vega export to multidimensional data ([4474f0c](https://github.com/biosimulations/biosimulations/commit/4474f0cf3cc7711b5e7b5dd9062e8d8e1c2585d9))
* expanded COMBINE archive creation to all types of model changes ([f65ffc5](https://github.com/biosimulations/biosimulations/commit/f65ffc5e6dde5a1edbbebbee83a20a73b1305479))


### Performance Improvements

* **api:** extract the combine archive directly on s3 via streams ([b5a0f08](https://github.com/biosimulations/biosimulations/commit/b5a0f0869610dc73e5119ecf44054e35f8354e15)), closes [#3094](https://github.com/biosimulations/biosimulations/issues/3094)


### Reverts

* "chore(deps): update dependency typescript to v4.5.4" ([f40932e](https://github.com/biosimulations/biosimulations/commit/f40932e6a41f56d94343ccad0af54d7368012572))
* "chore(deps): update typescript-eslint monorepo to v5.7.0" ([dee846f](https://github.com/biosimulations/biosimulations/commit/dee846f58f5989ef2f40733d97cbbff2b5f0722c))
* "chore(deps): update typescript-eslint monorepo to v5.7.0" ([16bcd80](https://github.com/biosimulations/biosimulations/commit/16bcd80c16b51ca2b3471b62f3fe10e0e2b448c7))


### BREAKING CHANGES

* **api:** the download project endpoint will fail for all previously submitted simulation
runs. Runs submitted prior to this change will not be retrievable by the api or applications

## [7.0.0](https://github.com/biosimulations/biosimulations/compare/v6.1.0...v7.0.0) (2021-12-02)


### Bug Fixes

* **api,dispatch-service,hsds:** fixed retrying in APIs ([b3109a1](https://github.com/biosimulations/biosimulations/commit/b3109a1fa5d37d0a65134adde6d9aaef90b19ac0))
* **config:** reverted changes to localhost ([04f2c42](https://github.com/biosimulations/biosimulations/commit/04f2c42263b8b98facda48a712bfb9f7e2d2e50d))
* corrected BioSimulators auth audience ([5c8b18d](https://github.com/biosimulations/biosimulations/commit/5c8b18d0d973bb39c8ba34375a9ba33219a76b44))
* corrected BioSimulators auth audience ([fce4347](https://github.com/biosimulations/biosimulations/commit/fce434727d36d28f048ea3b2d852ab56574a6768))
* corrected filtering for numerical columns ([d4bf35f](https://github.com/biosimulations/biosimulations/commit/d4bf35f09605fe5210c5b3ae953e0f6156624bec))
* **dispatch:** fix displaying of data in cases where simulation fails and metadata is not present ([da18704](https://github.com/biosimulations/biosimulations/commit/da18704a8b8fef837b9a835d12e95c6f6644d237)), closes [#3705](https://github.com/biosimulations/biosimulations/issues/3705)
* **dispatch:** fix file uploading ([30f9cb9](https://github.com/biosimulations/biosimulations/commit/30f9cb93e6ed9c52b08a62655de5693bd9552325)), closes [#3719](https://github.com/biosimulations/biosimulations/issues/3719)
* **dispatch:** fix uploading of files ([f597e90](https://github.com/biosimulations/biosimulations/commit/f597e9050b7385b8e17b7179476512d15cb3d723)), closes [#3719](https://github.com/biosimulations/biosimulations/issues/3719)
* fixed validation by upgrading from broken version of class-transformer ([ca20307](https://github.com/biosimulations/biosimulations/commit/ca203077941683837fb5426d1a7eb619d49d6b4c))
* **platform,ui:** fixed project browse for mobile ([fcfb552](https://github.com/biosimulations/biosimulations/commit/fcfb552a52b641d9f616f1e0785929d1b4bf413e))
* **platform:** fixed position of seach/filter button ([133f08c](https://github.com/biosimulations/biosimulations/commit/133f08c99c691b4cb62ba5c548abb2d7a41a7013))
* **ui:** fixed autocomplete filter ([cb26564](https://github.com/biosimulations/biosimulations/commit/cb26564c503e9d10f4423c21e27993e7c77e1802))


### Features

* added ability to filter and search projects ([21c14b4](https://github.com/biosimulations/biosimulations/commit/21c14b4fec4c1a0e56773b16e7b9e9e46c306636))
* added new component to enable components for routes to push buttons into the breadcrumbs area ([7918db1](https://github.com/biosimulations/biosimulations/commit/7918db1fccba56f3b9127bccdf042ddb7e7d83c4))
* **api,dispatch,dispatch-service:** expanded support for failed simulation runs ([9e7e71c](https://github.com/biosimulations/biosimulations/commit/9e7e71c80a0a81de3195648304456dfcf293b00c))
* **api,platform:** added owners, organizations to project view with hyperlinks ([63d9457](https://github.com/biosimulations/biosimulations/commit/63d9457a5244212523b727c3f7a01117b7253711))
* **api,platform:** started to display ownership of projects ([6f8378a](https://github.com/biosimulations/biosimulations/commit/6f8378a14d2f3d3291b36d11ada2581d1e6bd2be))
* **api:** add check for data service to status check ([d8fbbc5](https://github.com/biosimulations/biosimulations/commit/d8fbbc5308915d23e23dcd9b5962c869236f8988)), closes [#3649](https://github.com/biosimulations/biosimulations/issues/3649)
* **api:** added new scope for externally validating simulation runs ([e3bd698](https://github.com/biosimulations/biosimulations/commit/e3bd698c989c48cda49e242d2677859d928b9324))
* **api:** added URLs to accounts, organizations ([c6926b9](https://github.com/biosimulations/biosimulations/commit/c6926b95abbd8873bc4b56fbe77038b2cc2a76be))
* **api:** began to limit publication requests to model repositories ([e830002](https://github.com/biosimulations/biosimulations/commit/e8300028c8f8e22a2ef59938b7c83880167df5a7))
* **api:** working on restricting requests for publication with simulation run requests ([af7b6a0](https://github.com/biosimulations/biosimulations/commit/af7b6a0d46cedff4741d07502ccf8ef07fbbaa89))
* **dispatch:** clarified units of columns of simulation runs table ([3ae498e](https://github.com/biosimulations/biosimulations/commit/3ae498ecc1902c6fd688889427b0fdf28ad25b2e))
* improved table searching for data with accents ([d33617b](https://github.com/biosimulations/biosimulations/commit/d33617b73d1968f74cacaeae0bdd92e2a3aa3339))
* **platform:** added filter for publication status ([509eaa4](https://github.com/biosimulations/biosimulations/commit/509eaa4a431eb481cde071edc666ee7bdbe3961f))
* **platform:** scroll to top on opening projects search/filter ([26de2b6](https://github.com/biosimulations/biosimulations/commit/26de2b6f1470cc9fbb661447ac0542c2c7ea3483))
* **platform:** started to add filtering and searching for projects ([2505d7f](https://github.com/biosimulations/biosimulations/commit/2505d7f2bdca62c3b24d4b7252cb2513a9b5aa1b))
* **simulators:** expanded simulators filters ([1a76b7e](https://github.com/biosimulations/biosimulations/commit/1a76b7e5934c30793673ac8b2250c5b370ba2fa9))
* **ui,platform:** added autocomplete filter for attributes with many values ([049869f](https://github.com/biosimulations/biosimulations/commit/049869f919dcf9b714fc2d27389d1b1b6ce21971))
* **ui:** add custom caruousel component ([9a61d4f](https://github.com/biosimulations/biosimulations/commit/9a61d4f45592aa0b4c6d7205171a7b3a87b0330f))
* **ui:** replace npn-slider with custom component ([1cc79d5](https://github.com/biosimulations/biosimulations/commit/1cc79d54bc20114247ee6af36e2d3ade44c2f64d)), closes [#3706](https://github.com/biosimulations/biosimulations/issues/3706)


### Performance Improvements

* **api,dispatch,dispatch-service,simulators-api:** removed unnecessary return of new resources ([5fce07f](https://github.com/biosimulations/biosimulations/commit/5fce07f6b52385ab5d04c276b635f78aa06e2eea))
* **api:** don't return the log after creating ([5e37c6e](https://github.com/biosimulations/biosimulations/commit/5e37c6eca3aea9420467884f75f05558b6d30c0b)), closes [#3609](https://github.com/biosimulations/biosimulations/issues/3609)


### Reverts

* Revert "refactor(auth): removed auth/open endpoint" ([f983628](https://github.com/biosimulations/biosimulations/commit/f983628bd2b048b2a9e2b3875c4508f7a7b62746))
* "refactor: cleaned up building Angular apps" ([d21b2ed](https://github.com/biosimulations/biosimulations/commit/d21b2edc4530ff644d92ac6f71385dc83052fa3d))
* revert "refactor: organized endpoints configuration" ([a5e93c3](https://github.com/biosimulations/biosimulations/commit/a5e93c38a268995474525dc71e1fb72ef8fdf968)), closes [#3625](https://github.com/biosimulations/biosimulations/issues/3625)
* revert change to build front end apps ([5870fbf](https://github.com/biosimulations/biosimulations/commit/5870fbf953b95bd55d6def44aac7a5628c0a7265))


### BREAKING CHANGES

* **api:** The logs post endpoint no longer returns the log that was created. For that, use a
GET request after posting the log.

## [6.1.0](https://github.com/biosimulations/biosimulations/compare/v6.0.2...v6.1.0) (2021-11-14)


### Bug Fixes

* **auth:** fix import ([a7ba6b8](https://github.com/biosimulations/biosimulations/commit/a7ba6b87404b621d22c90c9fb37164b836d39f85))
* **auth:** handle case of no custom permissions ([1c3d760](https://github.com/biosimulations/biosimulations/commit/1c3d760f42ea23bdb3c26a42dce81b09189f7f92))
* corrected capitalization of BioSimulations ([3d981ee](https://github.com/biosimulations/biosimulations/commit/3d981ee737b852488dfd7ff8aba1317c11b7f236))
* debugged testing COMBINE API ([35672ce](https://github.com/biosimulations/biosimulations/commit/35672ceb0ce44474f8644c1bef71584208d778c7))
* debugged testing COMBINE API ([d18cb86](https://github.com/biosimulations/biosimulations/commit/d18cb864f29500229fa24666d9525bc228324476))
* debugged testing COMBINE API ([2a1e6b3](https://github.com/biosimulations/biosimulations/commit/2a1e6b3fc24fcdba88f7cc839d3655d6b4be80ce))
* debugged testing COMBINE API ([c7e1cb9](https://github.com/biosimulations/biosimulations/commit/c7e1cb9efc3daacaa6eddc1fd48806d25ec9927b))
* **dispatch:** corrected run URLs in check simulation run tool ([a1894fa](https://github.com/biosimulations/biosimulations/commit/a1894fae9b406a62fb4d10a30d01c9988d76a95a))
* **dispatch:** fixed dispatch simulation run view; closes [#3088](https://github.com/biosimulations/biosimulations/issues/3088) ([26a8d5a](https://github.com/biosimulations/biosimulations/commit/26a8d5a49b577c89eb04cd4b58314122f2d844c1))
* **dispatch:** fixed highlight.js import for log formatting ([0320cc2](https://github.com/biosimulations/biosimulations/commit/0320cc201838549a4d786c117571db132ddad600))
* fixed links, warnings ([ecb68fe](https://github.com/biosimulations/biosimulations/commit/ecb68febccdb4f17308aa62f0172adfcc7c68554))
* fixed python code highlighting ([ef88e12](https://github.com/biosimulations/biosimulations/commit/ef88e1248bbf2f9c8493c6f9a94a7c359d4497d0))
* fixed typos, added spelling exceptions ([7c77dc2](https://github.com/biosimulations/biosimulations/commit/7c77dc2db2b38563fea78118b1b65adca0510198))
* removed example with COMBINE archive that intentionally fails ([04b4c6e](https://github.com/biosimulations/biosimulations/commit/04b4c6ed1923ad2fa8ec1ae1f963e638cbe95807))


### Features

* **api,dispatch-service:** improved error messages and retrying ([a6b2693](https://github.com/biosimulations/biosimulations/commit/a6b26935a12e66a5bc3ab0b00713d8613c7a3f5d))
* **api:** improved reporting of errors with inconsistent data ([1c6c223](https://github.com/biosimulations/biosimulations/commit/1c6c2233c353f534efa0503f63b314adb9b97738))
* **dispatch-service:** add logging to processing posts to api ([1ca1900](https://github.com/biosimulations/biosimulations/commit/1ca19009cd33f6fbbb79ba9fb15779e09e20c1c8))
* **dispatch-service:** add retries for posting processing results ([fe9cddc](https://github.com/biosimulations/biosimulations/commit/fe9cddce4c846f15fdfd002bb8465e0e48bb3bd8)), closes [#3531](https://github.com/biosimulations/biosimulations/issues/3531)
* improved docs ([ab6722b](https://github.com/biosimulations/biosimulations/commit/ab6722b558a0b4e0d5aca265c84a1d0afd9f2558))
* improved docs ([209a421](https://github.com/biosimulations/biosimulations/commit/209a421c5c177600dd3e88c9469020db3c2fa51c))

## [6.0.2](https://github.com/biosimulations/biosimulations/compare/v6.0.1...v6.0.2) (2021-11-10)


### Bug Fixes

* **hsds:** update the hsds client ([419bfe9](https://github.com/biosimulations/biosimulations/commit/419bfe9a7ef4b81eea9400e2a4aa7587d734b4bf)), closes [#3317](https://github.com/biosimulations/biosimulations/issues/3317)

## [6.0.1](https://github.com/biosimulations/biosimulations/compare/v6.0.0...v6.0.1) (2021-11-10)


### Bug Fixes

* **api:** fixed IsImageDigest validator for non-strings ([53501b1](https://github.com/biosimulations/biosimulations/commit/53501b1eb3f323658aa23ae9007d625f212f6c6f))

## [6.0.0](https://github.com/biosimulations/biosimulations/compare/v5.9.0...v6.0.0) (2021-11-04)


### Bug Fixes

* **api:** build fix for new axios types ([e7ea984](https://github.com/biosimulations/biosimulations/commit/e7ea9849f72eefb1a051c316601e9463c3f74b1b))
* **dispatch-service:** add a temporary check for mistructured logs ([a67bfa1](https://github.com/biosimulations/biosimulations/commit/a67bfa1347b827cc2162080b0167412e33afa12f)), closes [#3482](https://github.com/biosimulations/biosimulations/issues/3482) [#3482](https://github.com/biosimulations/biosimulations/issues/3482)
* **dispatch-service:** remove ssl skip when downloading archive ([b4bd2c0](https://github.com/biosimulations/biosimulations/commit/b4bd2c0553bfa253002efb4063f3be30b321c15c)), closes [#3092](https://github.com/biosimulations/biosimulations/issues/3092)
* **dispatch,api,dispatch-service:** fixed data model for exceptions in simulation run logs ([191f1d3](https://github.com/biosimulations/biosimulations/commit/191f1d36ffc8824eb900b4571b937f0c9191b258))
* fix imports ([a076005](https://github.com/biosimulations/biosimulations/commit/a076005a50f89fa1f5cc76d66ededc44d72633ae))
* **simulators:** fixed text overflow of simulator test results ([7a8b11a](https://github.com/biosimulations/biosimulations/commit/7a8b11a85fd8a5975e39763cfb9247801ac8f5d6))


### Code Refactoring

* **api:** cleaned up simulation run files ([0213fb5](https://github.com/biosimulations/biosimulations/commit/0213fb508b66665b98c64c713ab146505695e2b9))


### Features

* added endpoints for getting summaries of projects ([2e9a54e](https://github.com/biosimulations/biosimulations/commit/2e9a54ebba20f9883bbf2c3dd6297db906da1285))
* **api:** added database model and validation for logs ([6a3a344](https://github.com/biosimulations/biosimulations/commit/6a3a344aa14fbd4952758f676160a9ba472b071e))
* **api:** added endpoints for getting individual SED elements; closes [#3439](https://github.com/biosimulations/biosimulations/issues/3439) ([8aa64fc](https://github.com/biosimulations/biosimulations/commit/8aa64fc205494d27c4f444c9d669031dd3ab49cc))
* **combine-api:** add dynamic module for combine api-client ([e7c2448](https://github.com/biosimulations/biosimulations/commit/e7c24486db6f7bf5b4cbc2b9e9c2c6a67a34bc1b)), closes [#3180](https://github.com/biosimulations/biosimulations/issues/3180)
* **datamodel:** add validation for image digests ([6106e54](https://github.com/biosimulations/biosimulations/commit/6106e54f5d54792b24fc5994a7ac1a42bed0c790))
* **dispatch-service:** enhanced tracking of processing results ([b4f01e3](https://github.com/biosimulations/biosimulations/commit/b4f01e3f4428b64ddb8aa3671fe738ff024515cc))
* **dispatch:** updated publication form, finished switching to Endpoints ([a5bef72](https://github.com/biosimulations/biosimulations/commit/a5bef7248174dec74154087b92bcc6005fa87726))
* **exceptions:** improve error handling ([4a3e8c7](https://github.com/biosimulations/biosimulations/commit/4a3e8c78dcb35756172da30ef803d2954f264bc6))
* **hsds:** handle transient hsds query failures ([f4a19f5](https://github.com/biosimulations/biosimulations/commit/f4a19f531c80ec50d324dec009132450dd74edb2)), closes [#3413](https://github.com/biosimulations/biosimulations/issues/3413)
* **ontology:** added parent/child relationships to ontology terms ([4107f1d](https://github.com/biosimulations/biosimulations/commit/4107f1d3d1f17a7509262caae998fd0222ebf0c3))
* **simulators-api:** add validation for api models ([ce2c5bb](https://github.com/biosimulations/biosimulations/commit/ce2c5bb0616600289e94dc16a7dcc17cfcb27dc4))
* **simulators,dispatch,platform:** added status bar to bottom of apps; closes [#3210](https://github.com/biosimulations/biosimulations/issues/3210) ([3630c23](https://github.com/biosimulations/biosimulations/commit/3630c23ef11325b4bb47012e4ca58ec7fefb6b7c))


### Reverts

* **dispatch-service:** revert using new config service to provide basepath ([e479d2e](https://github.com/biosimulations/biosimulations/commit/e479d2e531383e55ed64823b98b41d7be9130f7f))


### BREAKING CHANGES

* **api:** moves simulation run file information from 'Simulation Files' collection

## [5.9.0](https://github.com/biosimulations/biosimulations/compare/v5.8.0...v5.9.0) (2021-10-24)


### Bug Fixes

* **deps:** update to nx v13 ([48dbf7c](https://github.com/biosimulations/biosimulations/commit/48dbf7cfe4002aed9fcc06237f9bc995539573c2))
* **exceptions:** handle cases of payload too large errors ([522cec8](https://github.com/biosimulations/biosimulations/commit/522cec85bce117adbd24240ba026f705e782d185)), closes [nestjs/nest#5990](https://github.com/nestjs/nest/issues/5990) [#3349](https://github.com/biosimulations/biosimulations/issues/3349)
* **exceptions:** improve error handling ([083091e](https://github.com/biosimulations/biosimulations/commit/083091e8103340d6fd80534ec05cbc0d815e81cf))


### Features

* **api:** add caching to results and ontology endpoints ([bb2a991](https://github.com/biosimulations/biosimulations/commit/bb2a991c77ee32bb43ae9edb02938ca345700544))
* **api:** add caching to results endpoints ([ed54363](https://github.com/biosimulations/biosimulations/commit/ed5436370f4a9be2ed253d69417144d0db10e89b))
* **api:** add health check for job queue ([a05556e](https://github.com/biosimulations/biosimulations/commit/a05556ef33817253b32e52e8519806599e966723))
* **api:** add health module and endpoints ([da036f5](https://github.com/biosimulations/biosimulations/commit/da036f59d5f2479790a5d4827a2b830432ea35af))
* **api:** add various health checks and endpoints ([d30c7d2](https://github.com/biosimulations/biosimulations/commit/d30c7d25d62268f675cff1804b38309747c68e0b))
* **api:** setup results cache with REDIS ([3e40fde](https://github.com/biosimulations/biosimulations/commit/3e40fdebe38ca3e3765823bf7f98fd33895a7c9f))
* **dispatch-service:** add limit for retries of status for jobs ([f187b03](https://github.com/biosimulations/biosimulations/commit/f187b03024215b0cb2d53858085197b9d45736b9))
* **dispatch,platform:** added structured data for projects, simulation runs: ([7135ed0](https://github.com/biosimulations/biosimulations/commit/7135ed01f3107407906b81ad5b9ac94846f36f82))
* **dispatch,simulators:** added structured data tutorials ([5c23e8d](https://github.com/biosimulations/biosimulations/commit/5c23e8da1bcc2e12163f82aebf98c3a90c14bab8))
* **dispatch,simulators:** encoded FAQs into Schema.org ([936fb44](https://github.com/biosimulations/biosimulations/commit/936fb4499a2b132945b39fc2148428d26a268189))
* **exceptions:** dont process health check http exceptions ([b75747e](https://github.com/biosimulations/biosimulations/commit/b75747e6e04c511899037322e4829ae12433f3ab))
* **simulators-api,api:** added clearer payload too large messages ([98fb49d](https://github.com/biosimulations/biosimulations/commit/98fb49d651c01d621d175051cb030621b273034a))
* **simulators-api:** add health checks for simulators-api ([26d3b63](https://github.com/biosimulations/biosimulations/commit/26d3b63d37b1f271c0fc327b4a8e2c4981565651))
* **simulators:** added structured data for simulators as software applications ([7618db3](https://github.com/biosimulations/biosimulations/commit/7618db3ee539c31be2167a8d717cb4a5fac2c798))


### Reverts

* **simulators-api,api:** revert partially the changes in 98fb49d651c01d621d175051cb030621b273034a ([d88fcea](https://github.com/biosimulations/biosimulations/commit/d88fcead0eebe59f6394c9befca9e6ac7e132c83))

## [5.8.0](https://github.com/biosimulations/biosimulations/compare/v5.7.3...v5.8.0) (2021-10-20)


### Features

* **api:** enabling simulation run requests with latest version of a simulator ([2b7c2d9](https://github.com/biosimulations/biosimulations/commit/2b7c2d92ed10825bfce1a6c35a5a1b4908ceeebe))
* **config:** added endpoint for latest versions of simulators ([2a582ee](https://github.com/biosimulations/biosimulations/commit/2a582ee364a6308f42ae632b2487edd050964c87))
* **dispatch:** recorded simulator versions and digests for simulation runs ([f2abb39](https://github.com/biosimulations/biosimulations/commit/f2abb39790aa24c582c42521aca38f3dfbecaa56))
* **simulators:** added validation that version isn't reserved word 'latest' ([47843a8](https://github.com/biosimulations/biosimulations/commit/47843a8b687ff2ad001de286a07aa813735ca0c2))

## [5.7.3](https://github.com/biosimulations/biosimulations/compare/v5.7.2...v5.7.3) (2021-10-19)


### Bug Fixes

* **api:** correct field name to get values from dataservice ([53a6bbc](https://github.com/biosimulations/biosimulations/commit/53a6bbc5d2ccab9ea086f405656bf26d0cb2bacb)), closes [#3313](https://github.com/biosimulations/biosimulations/issues/3313)
* **api:** fix typo with checks ([caada0a](https://github.com/biosimulations/biosimulations/commit/caada0ab6e0f806508ee4d29966e5bcb0ca85109))

## [5.7.2](https://github.com/biosimulations/biosimulations/compare/v5.7.1...v5.7.2) (2021-10-19)


### Bug Fixes

* **api:** aligned parameter name in documentation ([0100a27](https://github.com/biosimulations/biosimulations/commit/0100a274725e22654ce4cc4534e5b0ff64a42de1))
* **simulators-api:** corrected put method; closes [#3305](https://github.com/biosimulations/biosimulations/issues/3305) ([57c34be](https://github.com/biosimulations/biosimulations/commit/57c34bef1f1370fdcd0a235020e4e7b5f20f5d54))

## [5.7.1](https://github.com/biosimulations/biosimulations/compare/v5.7.0...v5.7.1) (2021-10-18)


### Bug Fixes

* **account-api:** fixed route parameter names ([c7d8712](https://github.com/biosimulations/biosimulations/commit/c7d8712b0614d325bc1bf7fcc9145effcb6f2e60))
* **api:** fixed route parameter names ([1496f0f](https://github.com/biosimulations/biosimulations/commit/1496f0fc1c9ca0d63a620f73c216d9ab0752cce2))

## [5.7.0](https://github.com/biosimulations/biosimulations/compare/v5.6.2...v5.7.0) (2021-10-18)


### Bug Fixes

* **api:** fix docs and typing of open api definition ([#3307](https://github.com/biosimulations/biosimulations/issues/3307)) ([0640c6a](https://github.com/biosimulations/biosimulations/commit/0640c6aff0dd6a2d558d44be620d042b6e7ba49d)), closes [#3304](https://github.com/biosimulations/biosimulations/issues/3304)
* **api:** fix param name parsing for projectId ([af3c406](https://github.com/biosimulations/biosimulations/commit/af3c4069f352839a33696e3cf7c5dbb45d20c210))
* **combine-service:** added missing Swagger templates to Docker image ([f27b8c8](https://github.com/biosimulations/biosimulations/commit/f27b8c8831b11fffedd3355bc9668f78ad2e080c))


### Features

* **combine-service:** added health endpoint ([0d356d3](https://github.com/biosimulations/biosimulations/commit/0d356d347f11c3377f14ddbb66af8823198eaa30))
* **combine-service:** increased file upload limit, clarified error message ([5dfa25c](https://github.com/biosimulations/biosimulations/commit/5dfa25cca9fe73e2ef4b1af881f646bf9b022d5e))

## [5.6.2](https://github.com/biosimulations/biosimulations/compare/v5.6.1...v5.6.2) (2021-10-18)


### Bug Fixes

* **api:** fix module import ([d47c376](https://github.com/biosimulations/biosimulations/commit/d47c376915713fbc2ae97b296318061374b4bc10))
* **dispatch:** corrected file types for validate OMEX metadata form ([98794cd](https://github.com/biosimulations/biosimulations/commit/98794cd1397cb1ee921fad02b09d83ca82e7bd3f))

## [5.6.1](https://github.com/biosimulations/biosimulations/compare/v5.6.0...v5.6.1) (2021-10-18)


### Bug Fixes

* **api:** fix permissions for endpoints ([f00f6d1](https://github.com/biosimulations/biosimulations/commit/f00f6d116cc64525656f4b36030ee0109f9ff3b0)), closes [#3242](https://github.com/biosimulations/biosimulations/issues/3242)
* update client ids for api docs ([1ec36bb](https://github.com/biosimulations/biosimulations/commit/1ec36bb377ab939062da27455d199dbc3e4ada25))
## [5.6.0](https://github.com/biosimulations/biosimulations/compare/v5.5.0...v5.6.0) (2021-10-17)


### Bug Fixes

* **api:** update api to use updated hsds client ([e7832aa](https://github.com/biosimulations/biosimulations/commit/e7832aaf2c6d34c78903227dfc7e7db254b5b73d))
* **dispatch:** add flag to skip downloading test results to simulators service ([cbd63cc](https://github.com/biosimulations/biosimulations/commit/cbd63ccc6ecf757dec0da6cf000c200f2e8ebc4c)), closes [#3197](https://github.com/biosimulations/biosimulations/issues/3197)
* **dispatch:** fix alg list to empty list to prevent crashing ([338be99](https://github.com/biosimulations/biosimulations/commit/338be997320e053b054137a44546a82a21c13e34))
* make changes to update mongoose ([4653155](https://github.com/biosimulations/biosimulations/commit/46531556323e5d84cc96afbd87f4a58aa335be23))
* **platform,dispatch:** corrected link to simulation results in files tab ([00aa363](https://github.com/biosimulations/biosimulations/commit/00aa363a059fc9b0f0f95d316c1313efadb890b7))
* **platform:** redirect to 404 for non-existent projects; closes [#3234](https://github.com/biosimulations/biosimulations/issues/3234) ([173439a](https://github.com/biosimulations/biosimulations/commit/173439a2c322d933dc2ad4f5f1b3bcbeee666e80))
* **simulators:** fix json-ld metadata on index.html ([95b3d98](https://github.com/biosimulations/biosimulations/commit/95b3d983c1be2ddf633570214f8d35a920d98f1e))


### Features

* **api:** add custom styling to swagger ui ([90830c0](https://github.com/biosimulations/biosimulations/commit/90830c02f53df0fc9f86c7211b2fa3359a1b8807))
* **hsds:** update client library ([5746832](https://github.com/biosimulations/biosimulations/commit/574683246e7abb8a338d35d2a7a98ba25ce96d41))
* **simulators:** added repository digest to image model; closes [#3194](https://github.com/biosimulations/biosimulations/issues/3194) ([1293410](https://github.com/biosimulations/biosimulations/commit/129341077889888b41397cf471389d8921c8c2f8))
* **simulators:** expanded full text search; closes [#3209](https://github.com/biosimulations/biosimulations/issues/3209) ([c877189](https://github.com/biosimulations/biosimulations/commit/c877189fc5ece4943f3577602ff770660cdf01c0))


### Reverts

* **deps:** revert update dependency eslint to v8 ([0fdb3d8](https://github.com/biosimulations/biosimulations/commit/0fdb3d81a59710e4667974ecf0e42c8ec65ffd34))

## [5.5.0](https://github.com/biosimulations/biosimulations/compare/v5.4.0...v5.5.0) (2021-10-09)


### Bug Fixes

* **simulators-api:** add biosimulations.org to cors ([c2eea89](https://github.com/biosimulations/biosimulations/commit/c2eea89c81f266a00235bd1338102bbee7274dae))


### Features

* **api:** add case-insenstive unique index for project ids ([5f96f91](https://github.com/biosimulations/biosimulations/commit/5f96f91d83ae8835a313d6273eafab5e544f4e77)), closes [#3160](https://github.com/biosimulations/biosimulations/issues/3160)
* **api:** add controller level validation for project ids ([cdba9ce](https://github.com/biosimulations/biosimulations/commit/cdba9ce007fa92515f64147527824ca3df225808))
* **dispatch:** added check that simulation run was successful ([b4ade32](https://github.com/biosimulations/biosimulations/commit/b4ade32f2c3443d182f92ced5d12bf80341a260a))
* **platform:** added validation for project ids; closes [#3183](https://github.com/biosimulations/biosimulations/issues/3183) ([01b6178](https://github.com/biosimulations/biosimulations/commit/01b61788565e3c42c1b11acf1e2f02b349321cfc))

## [5.4.0](https://github.com/biosimulations/biosimulations/compare/v5.3.0...v5.4.0) (2021-10-08)


### Bug Fixes

* **combine-service:** update combine-service client to latest api changes ([74a70d8](https://github.com/biosimulations/biosimulations/commit/74a70d8fa73086277318b363dc58d1ebfa7d1970))
* **dispatch,platform:** fixed visualization rendering when groups of datasets are selected ([f205ab4](https://github.com/biosimulations/biosimulations/commit/f205ab457eab7712133d859e6ae4cfeb8dd16d63))
* **dispatch:** handle cases when metadata is empty without throwing error ([57a16c5](https://github.com/biosimulations/biosimulations/commit/57a16c54e496929a482fc69edc9bf1b2ca165e9d))


### Features

* **combine-service:** added endpoints for validation ([f602e65](https://github.com/biosimulations/biosimulations/commit/f602e651107f7efd195bb05ef41650bfe64009cd))
* **combine-service:** updated to biosimulators-utils 0.1.130 ([e46d13e](https://github.com/biosimulations/biosimulations/commit/e46d13ecc28e0e86ff2377068f8e11872276d327))
* **dispatch:** add some error handling ([6cc7f62](https://github.com/biosimulations/biosimulations/commit/6cc7f625451306c04e53884ab675d33ffd1fd5b8)), closes [#3088](https://github.com/biosimulations/biosimulations/issues/3088)
* **dispatch:** added forms for validating models, simulations and metadata ([d7991ed](https://github.com/biosimulations/biosimulations/commit/d7991ed6ce4d5cc14e35fbd70143ebb1d21705ca))
* **dispatch:** added options to project validation ([eebd16f](https://github.com/biosimulations/biosimulations/commit/eebd16fb26a66515952708e1f8e11638f234e93e))


### Reverts

* **combine-service:** reverted URL for COMBINE API ([d1dd8b5](https://github.com/biosimulations/biosimulations/commit/d1dd8b595ae1777c4e73b1e90d83f6a33603730c))

## [5.3.0](https://github.com/biosimulations/biosimulations/compare/v5.2.0...v5.3.0) (2021-10-07)


### Bug Fixes

* added 'master' attribute to file object ([5f4f722](https://github.com/biosimulations/biosimulations/commit/5f4f722d848e1e0a37802119b860b945677417a6))
* **api,mail-service:** update api client to return observable ([d26ce89](https://github.com/biosimulations/biosimulations/commit/d26ce8906de9275f40790ea4afbe6a97e404dd0d)), closes [#3102](https://github.com/biosimulations/biosimulations/issues/3102)
* **api:** add permissions to get all specs ([dbce421](https://github.com/biosimulations/biosimulations/commit/dbce4210625bc6ed6bb64eb7e3a170e385c88e55)), closes [#3136](https://github.com/biosimulations/biosimulations/issues/3136)
* **deps:** update dependency rxjs to v7.3.1 ([8cb2d32](https://github.com/biosimulations/biosimulations/commit/8cb2d3209fbcdc58771c4d5517aa5095376e87b5))
* **deps:** update nest ([f7a97e6](https://github.com/biosimulations/biosimulations/commit/f7a97e60127d9e66103dd331ce51038431398433))
* **dispatch:** fixed lint issue ([660bcb8](https://github.com/biosimulations/biosimulations/commit/660bcb8dcb6957bdcfc0d7469fcafb21f29f84c9))
* fixed lint issue ([e82e6fc](https://github.com/biosimulations/biosimulations/commit/e82e6fc06df38cc330dd7539e0cc7e18c45caab1))
* **platform,dispatch:** fixed plotly tests ([f2d8353](https://github.com/biosimulations/biosimulations/commit/f2d83535195d6d1716f2c3e9c3dc5d74aa310852))
* **simulators-api:** allow cors for biosimulatiors.dev ([1d93452](https://github.com/biosimulations/biosimulations/commit/1d93452b12174b444c5d53417521660dde72be56))
* **ui:** fixed display of errors with Vega visualizations ([1270955](https://github.com/biosimulations/biosimulations/commit/12709555792627db191d6794477f72ca6d81c7c4))
* **ui:** fixed Vega export for 1-d heatmaps ([cd5713d](https://github.com/biosimulations/biosimulations/commit/cd5713d73f14cddac1febaeb977ebb77a5fa0dff))


### Features

* **api:** create project endpoints ([d1b9fe7](https://github.com/biosimulations/biosimulations/commit/d1b9fe73c719358ace007c9d67a43c4a1d1c6810)), closes [#3067](https://github.com/biosimulations/biosimulations/issues/3067)
* **combine-service:** added options for validation of COMBINE archives ([42febbe](https://github.com/biosimulations/biosimulations/commit/42febbedc289d76f41bd1655014e1c9a172ff643))
* **combine-service:** added options to control COMBINE archive validation ([b4c0c12](https://github.com/biosimulations/biosimulations/commit/b4c0c123ea6ac776cc912e4eed964f773477d403))
* **combine-service:** added timeout for simulation execution ([8eb8deb](https://github.com/biosimulations/biosimulations/commit/8eb8debcab97b215260b9e6646751850b1a47593))
* **combine-service:** update combine-api client ([77c2f6d](https://github.com/biosimulations/biosimulations/commit/77c2f6df36cc59be74d4271938c34b5f074608a1))
* **datamodel:** add project datamodel ([57cf45c](https://github.com/biosimulations/biosimulations/commit/57cf45ceda0ba7d32f909e5d07fb7200e8dbbdee))
* **dispatch,platform,simulators:** improve recognition of Vega files by media type ([5f0051b](https://github.com/biosimulations/biosimulations/commit/5f0051bbfbd4e2ab6cd51d98a4352487ea94d640))
* **dispatch:** added options for validating COMBINE archives ([7d0a815](https://github.com/biosimulations/biosimulations/commit/7d0a8159798132fd65d93d7c91881a834b389c83))
* **platform:** added export of visualizations to Vega and COMBINE archives ([fce6731](https://github.com/biosimulations/biosimulations/commit/fce67312fa9432ed33ee0714a30d81d1496de111))
* **platform:** added heatmap, line plots ([84898db](https://github.com/biosimulations/biosimulations/commit/84898db3682a5b9321132a8a3058a6663c56113a))
* **platform:** added histogram visualization ([7a6abfa](https://github.com/biosimulations/biosimulations/commit/7a6abfafe6fad3dcebccff2977fec994839cdcbc))
* **platform:** added SED-ML visualizations ([2ac40b3](https://github.com/biosimulations/biosimulations/commit/2ac40b37cf40c713fd798c92fa761aa080586077))
* **platform:** added simulation types and algorithms to simulaton overview ([4c768e5](https://github.com/biosimulations/biosimulations/commit/4c768e594caee33d7cee79f157eaa70ba6703a3c))
* **platform:** added vega export for heatmaps and line plots ([a8c9ad4](https://github.com/biosimulations/biosimulations/commit/a8c9ad418407d7ec7d7f42e21ddf6a0d0dd5d966))
* **platform:** front end displays projects from api ([9ecfa80](https://github.com/biosimulations/biosimulations/commit/9ecfa80ce7f70a45762563b88c82c0bf0e0cf3a0)), closes [#3149](https://github.com/biosimulations/biosimulations/issues/3149)
* **ui:** added ability to attach hyperlinks to menu items ([7ca3d10](https://github.com/biosimulations/biosimulations/commit/7ca3d10209fa04d77cbc220664f0d1870c542c12))


### Reverts

* **deps:** revert 235c9db3e9649cdb8b42e6575517aa651f9e1c2d ([05cb6f3](https://github.com/biosimulations/biosimulations/commit/05cb6f39142d576e4d9c866a6ee2177589aadbd3))

## [5.2.0](https://github.com/biosimulations/biosimulations/compare/v5.1.1...v5.2.0) (2021-10-04)


### Bug Fixes

* **api:** add authentication to post metadata ([999e2b9](https://github.com/biosimulations/biosimulations/commit/999e2b9e11ddaa26e001d89ccbbeb26541cc3f5a)), closes [#2865](https://github.com/biosimulations/biosimulations/issues/2865)
* **deps:** update dependency @sendgrid/mail to v7.4.7 ([bae8e96](https://github.com/biosimulations/biosimulations/commit/bae8e968c0062da5f477ca514e3ed3c08e6686e3))
* **dispatch-service:** fix processing of environment variables ([cb0aa04](https://github.com/biosimulations/biosimulations/commit/cb0aa04cbdc7bd8f2a3e7cf818c541470d7c4519))
* **dispatch-service:** process metadata with other processing ([7f8f44a](https://github.com/biosimulations/biosimulations/commit/7f8f44aecc741dcd746e0de71489bcd5a9c68319)), closes [#3046](https://github.com/biosimulations/biosimulations/issues/3046)


### Features

* **api:** add endpoints to get particular specifications for simulation runs ([87eede1](https://github.com/biosimulations/biosimulations/commit/87eede1c0283bdc3a4a5b6614dba06354a96ddf6))
* **api:** add file object ([39f25f3](https://github.com/biosimulations/biosimulations/commit/39f25f3668d7e249ae1b9834b6ef66cb58350208)), closes [#2914](https://github.com/biosimulations/biosimulations/issues/2914)
* **api:** create specifications object and endpoints ([aca4786](https://github.com/biosimulations/biosimulations/commit/aca4786282526813eab3280c87d495617c7a2ef8))
* **dispatch-service:** process files and sedml specs ([fb37624](https://github.com/biosimulations/biosimulations/commit/fb376243c33043d4945a2849d0df50a54332fcda))
* **dispatch-service:** send sedml specifications to the api ([7d6a8c7](https://github.com/biosimulations/biosimulations/commit/7d6a8c7baec360f985459534f9fd5d67e4342260))

## [5.1.1](https://github.com/biosimulations/biosimulations/compare/v5.1.0...v5.1.1) (2021-10-01)


### Bug Fixes

* **platform:** corrected handling of software license keys for simulation ([84bef20](https://github.com/biosimulations/biosimulations/commit/84bef20d9c035168d7d74fcffdeafae97091af8f))
* **platform:** corrected handling of software license keys for simulation ([aba13fe](https://github.com/biosimulations/biosimulations/commit/aba13fe9c30c2355d6d63a9a3546fa78062b82d6))

## [5.1.0](https://github.com/biosimulations/biosimulations/compare/v5.0.0...v5.1.0) (2021-09-29)


### Bug Fixes

* **deps:** update dependency auth0 to v2.36.2 ([#3076](https://github.com/biosimulations/biosimulations/issues/3076)) [skip ci] ([a696fbc](https://github.com/biosimulations/biosimulations/commit/a696fbcc9f495fef5b20fba7235a8fb7d590b4d9))
* **dispatch-api:** fix cors for biosimulators ([22439d2](https://github.com/biosimulations/biosimulations/commit/22439d2d6d2c16ce71dc652e5353a9f902a29c6a))


### Features

* **dispatch:** require configuration of academic use for commercial solvers ([6c5307c](https://github.com/biosimulations/biosimulations/commit/6c5307c31dad3f59b2e7c01b07ea1b83218a7ff0))

## [5.0.0](https://github.com/biosimulations/biosimulations/compare/v4.6.0...v5.0.0) (2021-09-28)


### Bug Fixes

* **dispatch-service:** fix type error when processing metadata ([a646c98](https://github.com/biosimulations/biosimulations/commit/a646c98ef2433b5b20b754fd637f326513aa57e1))
* make datamodel consistent for license ([4b95e4d](https://github.com/biosimulations/biosimulations/commit/4b95e4d89af83b69334602f53ffacaa0744e5aff)), closes [#3050](https://github.com/biosimulations/biosimulations/issues/3050)
* **dispatch-api:** remove extra slash for metadata uris ([b627e74](https://github.com/biosimulations/biosimulations/commit/b627e7491995ae1c1e70feea93b1e7f4cc53902a)), closes [#3052](https://github.com/biosimulations/biosimulations/issues/3052)
* ensure external url is used for combine api ([2d98aba](https://github.com/biosimulations/biosimulations/commit/2d98aba10be03163e9b43aa67da69a95910eb763))
* **dispatch:** corrected when metadata about simulation projects is retrieved ([1682d83](https://github.com/biosimulations/biosimulations/commit/1682d83d6870e4bcfa1d4b895d3a88d2cee60285))


### Code Refactoring

* consolidate backend apis ([b27bd0e](https://github.com/biosimulations/biosimulations/commit/b27bd0e260336df3553b1b3a7e3447c0e26ac716)), closes [#2724](https://github.com/biosimulations/biosimulations/issues/2724)


### Features

* **dispatch-api:** ensure only public models are shown for platform ([#3045](https://github.com/biosimulations/biosimulations/issues/3045)) ([5619c03](https://github.com/biosimulations/biosimulations/commit/5619c03d7fc088d0b3be33136935c80e7cb9c862)), closes [#3044](https://github.com/biosimulations/biosimulations/issues/3044)
* **dispatch-api:** extract files to s3 and replace combine archive file extraction endpoint ([56f8413](https://github.com/biosimulations/biosimulations/commit/56f84133b193c4d54f77c33ba2c01105df6162e3)), closes [#2945](https://github.com/biosimulations/biosimulations/issues/2945)
* **dispatch-api:** upload omex files to s3 from url ([4d8f780](https://github.com/biosimulations/biosimulations/commit/4d8f78058d3ea7050d913ad651c907c0a631a3f4))


### Reverts

* **dispatch-api:** revert permissions change in 3175c6378160f34e8389b6e501ea2534eb9d4c12 ([5ab7d08](https://github.com/biosimulations/biosimulations/commit/5ab7d08cc922778436646dca0331aad7bafef0d3))


### BREAKING CHANGES

* The ontology, dispatch, and platform apis are consolidated into one main backend
api for biosimulations. There is a seperate api for biosimulators. The combine-service also provides a rest api that is mostly intended for internal use.

## [4.6.0](https://github.com/biosimulations/biosimulations/compare/v4.5.0...v4.6.0) (2021-09-27)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.995.0 ([19b509a](https://github.com/biosimulations/biosimulations/commit/19b509a7100e59dbcdf3c2674262ec2bf4333f77))
* **dispatch:** improved handling of undefined simulations in simulations browse view; closes [#2999](https://github.com/biosimulations/biosimulations/issues/2999) ([bf47994](https://github.com/biosimulations/biosimulations/commit/bf4799404320d07437797804973e933d7f147c7e))
* **platform:** handle cases when metadata is missing ([f738a74](https://github.com/biosimulations/biosimulations/commit/f738a741caba0da783e522c8d74e86c04b9e4aa8))


### Features

* **simulators:** improved specification of simulator CLIs; closes [#3015](https://github.com/biosimulations/biosimulations/issues/3015) ([c396bfa](https://github.com/biosimulations/biosimulations/commit/c396bfa79a7f9abbc0a3e6b374f96d726fe5eaa7))
* create alternate vega view component ([453a109](https://github.com/biosimulations/biosimulations/commit/453a1094677a18ccc2019bd953d618e4e81c6e2f))
* **dispatch:** added license confirmation to publish run form ([1889e1e](https://github.com/biosimulations/biosimulations/commit/1889e1ec4f50579c51838494f7d693fdb5a87415))
* **dispatch,platform:** updated terms about granting BioSimulations to distribute projects ([6ee2bb1](https://github.com/biosimulations/biosimulations/commit/6ee2bb1779e3f02df9dee82bde75e11fa5dc5f9d))
* **platform:** add code to get metadata and specs of visualiazations ([0cdc472](https://github.com/biosimulations/biosimulations/commit/0cdc4727b201858adec4cb40565bb643c01c2905))
* **platform:** add support for showing vega figures ([f011480](https://github.com/biosimulations/biosimulations/commit/f011480096ee29f42e8aab833dbaec89a1bc1596))
* **ui:** allow for conditional loading of tabs ([1b4779f](https://github.com/biosimulations/biosimulations/commit/1b4779f177748922805e8aaad196d0ef412a5239))
* updated biosimulators-utils, biosimulators-bionetgen ([ce75ea6](https://github.com/biosimulations/biosimulations/commit/ce75ea6a46e880b1c021a336992d7acc311f200d))


### Performance Improvements

* **platform:** fix tests ([a8f6d0e](https://github.com/biosimulations/biosimulations/commit/a8f6d0ed968a9ccccd45249db4b71313ef7a5f5b))

## [4.5.0](https://github.com/biosimulations/biosimulations/compare/v4.4.2...v4.5.0) (2021-09-26)


### Features

* **simulators:** improved simulator usage examples; closes [#3016](https://github.com/biosimulations/biosimulations/issues/3016) ([a56dea2](https://github.com/biosimulations/biosimulations/commit/a56dea24e33e33b7db215a0cc1f62f4bc94dac7d))

## [4.4.2](https://github.com/biosimulations/biosimulations/compare/v4.4.1...v4.4.2) (2021-09-23)


### Bug Fixes

* **simulators-api:** fixed sorting of simulator versions with > 4 points ([020073a](https://github.com/biosimulations/biosimulations/commit/020073ab0a9e8eb0344d7ab4ffee287a323caaed)), closes [#3008](https://github.com/biosimulations/biosimulations/issues/3008)


### Reverts

* **simulators-api:** revert [#3009](https://github.com/biosimulations/biosimulations/issues/3009) to prevent container crashing ([c7ff2cc](https://github.com/biosimulations/biosimulations/commit/c7ff2cc018e28dad55f771c6e778193082e9c2ff)), closes [#3008](https://github.com/biosimulations/biosimulations/issues/3008)

## [4.4.1](https://github.com/biosimulations/biosimulations/compare/v4.4.0...v4.4.1) (2021-09-22)


### Bug Fixes

* fixed sorting of simulator versions with > 4 points ([e0b60ce](https://github.com/biosimulations/biosimulations/commit/e0b60ce7e857a949d73f4daed4cf51e9e2b0ca91))

## [4.4.0](https://github.com/biosimulations/biosimulations/compare/v4.3.0...v4.4.0) (2021-09-22)


### Features

* **ontology:** updated to KiSAO 2.29 ([31e7d7b](https://github.com/biosimulations/biosimulations/commit/31e7d7ba7115c9f381b636113221edca9010397b))
* **platform:** improve styling of platform browse ([28e8e1d](https://github.com/biosimulations/biosimulations/commit/28e8e1d6858363325ae29ae3d71e1dea2a1b19c9))


### Reverts

* remove ui commit scope [skip ci] ([5051250](https://github.com/biosimulations/biosimulations/commit/50512504bd98cf55f2111ad0bf074ae5837260ea))

## [4.3.0](https://github.com/biosimulations/biosimulations/compare/v4.2.0...v4.3.0) (2021-09-15)


### Bug Fixes

* **deps:** update dependency @stoplight/json-ref-resolver to v3.1.3 ([#2986](https://github.com/biosimulations/biosimulations/issues/2986)) ([fe9c5f3](https://github.com/biosimulations/biosimulations/commit/fe9c5f372c69c2f8a11d0d401a312ecdb7bf3338))
* **deps:** update dependency bull to v3.29.2 ([#2987](https://github.com/biosimulations/biosimulations/issues/2987)) ([98cbebf](https://github.com/biosimulations/biosimulations/commit/98cbebf7752e6fd6b034e0a459b6abc9de106247))
* **dispatch:** proceed if metadata is missing ([#2998](https://github.com/biosimulations/biosimulations/issues/2998)) ([f09c633](https://github.com/biosimulations/biosimulations/commit/f09c633a4b69b096c5bd07a1377c451ea3ae3aa1)), closes [#2994](https://github.com/biosimulations/biosimulations/issues/2994)


### Features

* **simulators:** expanded specs for simulators ([32b100b](https://github.com/biosimulations/biosimulations/commit/32b100be3f5856bc8131427155031dc5abe1013a))
* expanded simulator specs ([f281cd3](https://github.com/biosimulations/biosimulations/commit/f281cd31b7c3c7c08dfc47162944dcfdbb7c4761))

## [4.2.0](https://github.com/biosimulations/biosimulations/compare/v4.1.0...v4.2.0) (2021-09-11)


### Bug Fixes

* **deps:** update dependency axios to v0.21.4 ([#2952](https://github.com/biosimulations/biosimulations/issues/2952)) ([c24c0e2](https://github.com/biosimulations/biosimulations/commit/c24c0e2eb96d6ff5c1a77f9408dee9c80e02c0f1))
* **deps:** update dependency ssh2 to v1.4.0 ([#2953](https://github.com/biosimulations/biosimulations/issues/2953)) ([a3afb2b](https://github.com/biosimulations/biosimulations/commit/a3afb2b7c862b7cf4bb86a451dce7380a64afa39))
* specify global setTimeout instead of window ([432b90a](https://github.com/biosimulations/biosimulations/commit/432b90a00f34c604bc5f2fd7038a925c52ea4fec))
* **dispatch-service:** restore check for empty env variables ([3f714ba](https://github.com/biosimulations/biosimulations/commit/3f714baf82ab04673fb33a28cc9f7daa9899c39b))
* restore mkdocs file location ([41a269e](https://github.com/biosimulations/biosimulations/commit/41a269e37f9c3e704baa614b4ab8de30d8ed6546))


### Features

* **account-api:** replace typegoose with mongoose ([69911a3](https://github.com/biosimulations/biosimulations/commit/69911a35ab0ca6c324ba53a097c814d53730d491))
* **dispatch-api:** handle errors and timeouts on uploa> ([dad35ea](https://github.com/biosimulations/biosimulations/commit/dad35ea5718abac450c682970e332ee4292890f3)), closes [#2860](https://github.com/biosimulations/biosimulations/issues/2860)
* **dispatch-service:** added passing software licenses from deployment secrets to Singularity run ([cc19999](https://github.com/biosimulations/biosimulations/commit/cc199990c6255b6b70bc436a9d16051602c4d0c5))
* **storage:** add simulation storage service and timeout for s3 uploads ([0c24173](https://github.com/biosimulations/biosimulations/commit/0c241737989289f79b13ebc8efa265cc7c6fc91f))

## [4.1.0](https://github.com/biosimulations/biosimulations/compare/v4.0.1...v4.1.0) (2021-09-06)


### Features

* **simulators:** added attribute to track installation instructions for Python APIs ([cb2b415](https://github.com/biosimulations/biosimulations/commit/cb2b415b4513170bbb09140e6cd9bff4b970d3ed))

## [4.0.1](https://github.com/biosimulations/biosimulations/compare/v4.0.0...v4.0.1) (2021-09-06)


### Bug Fixes

* **dispatch:** simulation results URLs for data visualizations ([9b7b879](https://github.com/biosimulations/biosimulations/commit/9b7b8795ac2b524b9089fad5a2c7916e3d1214f4))


### Reverts

* 39a60b17d640b62639f6594024f4ba4c66baedc5 ([f804cce](https://github.com/biosimulations/biosimulations/commit/f804cce9e3b3787a11b2989743e86407a4c014dd)), closes [#2959](https://github.com/biosimulations/biosimulations/issues/2959)

## [4.0.0](https://github.com/biosimulations/biosimulations/compare/v3.20.0...v4.0.0) (2021-09-04)


### Bug Fixes

* **dispatch:** properly encode uri to allow for fetching results ([dcbf044](https://github.com/biosimulations/biosimulations/commit/dcbf04433f7e105a01f501f3aa7172c82807ea41))


### Features

* update example simulation runs ([395f513](https://github.com/biosimulations/biosimulations/commit/395f513657c662d2b26b3d3b0de95cdd860ea326)), closes [#2951](https://github.com/biosimulations/biosimulations/issues/2951)


### BREAKING CHANGES

* simulation runs sumbitted prior to the update will not display on the dispatch app

## [3.20.0](https://github.com/biosimulations/biosimulations/compare/v3.19.0...v3.20.0) (2021-09-04)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.983.0 ([#2947](https://github.com/biosimulations/biosimulations/issues/2947)) ([39a60b1](https://github.com/biosimulations/biosimulations/commit/39a60b17d640b62639f6594024f4ba4c66baedc5))
* **dispatch-service:** correct the determination of the environment ([ce46d3b](https://github.com/biosimulations/biosimulations/commit/ce46d3bddb05195dd29408d05c00de53186336ae))
* fix default environment to dev ([c970ccd](https://github.com/biosimulations/biosimulations/commit/c970ccde7a533bb0db3f0ec8334308d3a1ee237d))
* new endpoint implementation ([ed42b6b](https://github.com/biosimulations/biosimulations/commit/ed42b6b27fbba4a97b8a931d6d9771b1563eddb9)), closes [#2943](https://github.com/biosimulations/biosimulations/issues/2943) [#2861](https://github.com/biosimulations/biosimulations/issues/2861) [#2859](https://github.com/biosimulations/biosimulations/issues/2859)


### Features

* **dispatch,dispatch-api:** move thumbnail processing to backend ([4495d6d](https://github.com/biosimulations/biosimulations/commit/4495d6d70e8fcb168fd5b4a38f70850171908d7b))
* **platform:** add page to view projects on platform ([e568d0c](https://github.com/biosimulations/biosimulations/commit/e568d0c16a5e67198e921d8705e45f137c680df2))


### Reverts

* revert commit 5dad745d1df0ffc3fb2fba8fc3b99b21b69b0521 ([f8cdd5b](https://github.com/biosimulations/biosimulations/commit/f8cdd5b338ec3fd7f8b7faf607ce893a9d343075))

## [3.19.0](https://github.com/biosimulations/biosimulations/compare/v3.18.0...v3.19.0) (2021-09-02)


### Features

* **combine-service:** updated to biosimulators-utils 0.1.115, biosimulators-amici 0.1.18 ([9ad2945](https://github.com/biosimulations/biosimulations/commit/9ad29450a51b8ff181a00fe57c70b660dc917a60))

## [3.18.0](https://github.com/biosimulations/biosimulations/compare/v3.17.0...v3.18.0) (2021-09-01)


### Bug Fixes

* **deps:** update dependency form-data to v4 ([#2925](https://github.com/biosimulations/biosimulations/issues/2925)) ([36a79ba](https://github.com/biosimulations/biosimulations/commit/36a79baccabb13e181e9ca5ee0cd2d0bff629697))
* **deps:** update dependency jwks-rsa to v2 ([#2928](https://github.com/biosimulations/biosimulations/issues/2928)) ([f4a3f10](https://github.com/biosimulations/biosimulations/commit/f4a3f107f5e43831e4cf72ffc63fcaf67a0026e3))
* **deps:** update dependency ssh2 to v1 ([#2930](https://github.com/biosimulations/biosimulations/issues/2930)) ([11111c7](https://github.com/biosimulations/biosimulations/commit/11111c76a5c943741397d3110189ac0d5ee53a86))


### Features

* **combine-service:** fixed error handling for run sim, simplified run sim options ([5e63d49](https://github.com/biosimulations/biosimulations/commit/5e63d49eb5a1ad5c27ac09dc970093f04ff79980))
* **dispatch:** added support for new SBO modeling framework terms ([80ee759](https://github.com/biosimulations/biosimulations/commit/80ee759d6be92545b01f999c1a7c0630fa43f43d))

## [3.17.0](https://github.com/biosimulations/Biosimulations/compare/v3.16.0...v3.17.0) (2021-09-01)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.980.0 ([#2906](https://github.com/biosimulations/Biosimulations/issues/2906)) ([163191d](https://github.com/biosimulations/Biosimulations/commit/163191d5cc1e24dbeb0440681761e175b633a759))


### Features

* add shared config file support ([976e578](https://github.com/biosimulations/Biosimulations/commit/976e57846a8c43fa10f8be4e70a8a1989bde683c))
* **dispatch:** call the metadata endpoint to get simulation metadata ([ae1054f](https://github.com/biosimulations/Biosimulations/commit/ae1054f6f101b170cb1408d13ffdcbb39f0b25a1)), closes [#2866](https://github.com/biosimulations/Biosimulations/issues/2866)

## [3.16.0](https://github.com/biosimulations/Biosimulations/compare/v3.15.0...v3.16.0) (2021-08-31)


### Bug Fixes

* **deps:** update dependency class-validator to v0.13.1 ([#2894](https://github.com/biosimulations/Biosimulations/issues/2894)) ([3676e59](https://github.com/biosimulations/Biosimulations/commit/3676e59d26a91d81bd7b12cd0287d619fa8e89ab))
* **deps:** update dependency nats to v2.2.0 ([#2886](https://github.com/biosimulations/Biosimulations/issues/2886)) ([7097edc](https://github.com/biosimulations/Biosimulations/commit/7097edc10eef06edc70ede51b1c2c4dc9eec5810))
* **deps:** update dependency rxjs to v7.3.0 ([#2895](https://github.com/biosimulations/Biosimulations/issues/2895)) ([c01604b](https://github.com/biosimulations/Biosimulations/commit/c01604bfcf4df5516158b442435f63002ef9720c))
* **deps:** update dependency stackdriver-errors-js to v0.10.0 ([43451aa](https://github.com/biosimulations/Biosimulations/commit/43451aa19accdc7c670a17cd6fec1aa660934234))
* **deps:** update dependency tslib to v2.3.1 ([#2888](https://github.com/biosimulations/Biosimulations/issues/2888)) ([fc78756](https://github.com/biosimulations/Biosimulations/commit/fc787565990a70489c6ee0a6b9450af7e92eb118))


### Features

* **dispatch:** added dry run option to example simulation submission ([1487880](https://github.com/biosimulations/Biosimulations/commit/1487880a87f82f48e29dafcfe9741bf9ff862cb7))
* **dispatch:** added example simulation run for RBApy ([277772e](https://github.com/biosimulations/Biosimulations/commit/277772e9d062f4a82fab455efbc2c67d088770ad))
* **dispatch-api:** set uris for metadata elements ([96e94fe](https://github.com/biosimulations/Biosimulations/commit/96e94fef4ea1a25ad95cba100c9b635618d29e7e))
* **simulators:** added ability to capture Python APIs in simulator specs ([5ed44cb](https://github.com/biosimulations/Biosimulations/commit/5ed44cb904349e4f91eb6b8cf62f264eb6eeebdd))

## [3.15.0](https://github.com/biosimulations/Biosimulations/compare/v3.14.0...v3.15.0) (2021-08-29)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.978.0 ([#2883](https://github.com/biosimulations/Biosimulations/issues/2883)) ([ca487aa](https://github.com/biosimulations/Biosimulations/commit/ca487aab71671f9a4a71668e253800fdc7e98708))
* **deps:** update dependency bull to v3.29.1 ([#2884](https://github.com/biosimulations/Biosimulations/issues/2884)) ([e6589e7](https://github.com/biosimulations/Biosimulations/commit/e6589e7d70baa6c9d6797d4fde95687d4798c19b))
* **deps:** update nest ([#2835](https://github.com/biosimulations/Biosimulations/issues/2835)) ([0e65b60](https://github.com/biosimulations/Biosimulations/commit/0e65b6084e75aae31ff8d089e6321f581fd8742d))


### Features

* **combine-service:** updated to Biosimulators-utils with support for RBA models ([610225b](https://github.com/biosimulations/Biosimulations/commit/610225b46f8bd9ed26c1ca632f01c051d5765dc8))
* updated biosimulators documentation links to docs.biosimulatos.org ([bfa49bb](https://github.com/biosimulations/Biosimulations/commit/bfa49bb7530cb656689eb8632b365110fb6b5aca))
* updated SBO for term for RBA ([fc64418](https://github.com/biosimulations/Biosimulations/commit/fc64418993cb06a162884461646b586801b7f37e))

## [3.14.0](https://github.com/biosimulations/Biosimulations/compare/v3.13.0...v3.14.0) (2021-08-27)


### Bug Fixes

* **combine-service:** dont change field "abstract" to "_abstract" ([591b0db](https://github.com/biosimulations/Biosimulations/commit/591b0db5bef0c4a260fb6f0583ed82b0eaf77481))
* **combine-service:** fix api client implementation ([d313f1c](https://github.com/biosimulations/Biosimulations/commit/d313f1cae8f64d1e93a0741bab2545d9b14d1109))
* **dispatch-service:** dont log error if job is not yet present ([b80db71](https://github.com/biosimulations/Biosimulations/commit/b80db71ed633cbf9c6cb2103a513d7dd1fc397cd))
* **dispatch-service:** fix url for posting metadata ([a5926e7](https://github.com/biosimulations/Biosimulations/commit/a5926e70d24985d0d50adfc65acdd4d6bbb8ca7d))
* **platform:** correct url for metadata ([c2a7a63](https://github.com/biosimulations/Biosimulations/commit/c2a7a63f57072598f9890905697edcfb461742f4))
* **platform:** fix unterminated string literal ([fe0343d](https://github.com/biosimulations/Biosimulations/commit/fe0343dc6fe808023c822b3ffeb004c094be96a2))


### Features

* **combine-service:** isolating simulation execution into separate processes ([9d23a5d](https://github.com/biosimulations/Biosimulations/commit/9d23a5d12571708ecf7b44bc83675cb5b5802a98))
* **combine-service:** update combine api client ([c1bb566](https://github.com/biosimulations/Biosimulations/commit/c1bb56636207d966bbd41f2dcf697008c43a27cf))
* **dispatch-service:** create handler to extract metadata ([16b469a](https://github.com/biosimulations/Biosimulations/commit/16b469a479061d643ffa8541a5db55e95f687404))
* **dispatch-service:** process and create metadata for simulation runs ([99df19d](https://github.com/biosimulations/Biosimulations/commit/99df19d50c786992aad6d35da1443fa3b8126c98))
* **exceptions:** change database errors to return 500 errors instead of 400 ([7390b7f](https://github.com/biosimulations/Biosimulations/commit/7390b7f3e35c6c1f0d1dfb8f145cf4ad8e5545c9))

## [3.13.0](https://github.com/biosimulations/Biosimulations/compare/v3.12.0...v3.13.0) (2021-08-24)


### Features

* **ontology:** updated to KiSAO 2.27 ([09ede72](https://github.com/biosimulations/Biosimulations/commit/09ede7208ce04542e5243ec96b4480471e98f8eb))

## [3.12.0](https://github.com/biosimulations/Biosimulations/compare/v3.11.0...v3.12.0) (2021-08-23)


### Bug Fixes

* **dispatch-api:** fields paramter is optional ([c3863e3](https://github.com/biosimulations/Biosimulations/commit/c3863e3c92709fb094ea6a9712aacbe59cdd412b))


### Features

* **dispatch-api:** implement metadata endpoints ([9d067e9](https://github.com/biosimulations/Biosimulations/commit/9d067e983cd625a8d706bc1cb3cfa2033bdabf62))

## [3.11.0](https://github.com/biosimulations/Biosimulations/compare/v3.10.0...v3.11.0) (2021-08-22)


### Features

* **combine-service:** enabled NEURON, NetPyNe for simulation ([19661df](https://github.com/biosimulations/Biosimulations/commit/19661df2b96697934d2f9ca3b9949cac1570554e))

## [3.10.0](https://github.com/biosimulations/Biosimulations/compare/v3.9.0...v3.10.0) (2021-08-21)


### Features

* **combine-service:** added endpoint for low latency simulation ([b44f5e3](https://github.com/biosimulations/Biosimulations/commit/b44f5e30df6b0e23b40ed302fa6a7a61a6969c09))
* **combine-service:** added endpoint for low-latency simulations ([931f3c7](https://github.com/biosimulations/Biosimulations/commit/931f3c7d25b329f55578cb0e7ba5e59ce7c91858))
* **combine-service:** added options to export simulation results in HDF5, zip formats ([c0cd699](https://github.com/biosimulations/Biosimulations/commit/c0cd6998dc850c5afbcee00fa2afec7390ffd235))
* **combine-service:** added simulator name to get simulaton tools endpoint ([1f708dd](https://github.com/biosimulations/Biosimulations/commit/1f708dd70e138a2c7d662b8689ef73566180e3ba))
* **combine-service:** added test to verify simulator APIs ([0bcaff3](https://github.com/biosimulations/Biosimulations/commit/0bcaff3f11fbc9e85285ef709432219d80fb1ef1))
* **combine-service:** pre-compiled Python code for faster initial calls ([f368de9](https://github.com/biosimulations/Biosimulations/commit/f368de9d973274c0e218e24075665a255d698e7f))
* **combine-service:** updated AMICI, GillesPy2, libSBMLSim ([fe8df2a](https://github.com/biosimulations/Biosimulations/commit/fe8df2a7b1da8c955210be78432565f60080d5bd))
* **combine-service:** updated to KiSAO 2.26, BioSimalators-utils 0.1.105 ([05b9dc5](https://github.com/biosimulations/Biosimulations/commit/05b9dc5ebe393c1e4c4143a0125867a3a79a58ea))
* **combine-service:** updated to KiSAO 2.26, BioSimulators-utils 0.1.105 ([74d0b36](https://github.com/biosimulations/Biosimulations/commit/74d0b369a933ec4b3f7b46c4028addfb794ef90d))
* **dispatch:** added simulator names to simulation tools menu in run form ([5170810](https://github.com/biosimulations/Biosimulations/commit/51708108ce772851690a4a3edff556513feb9368))
* **ontology:** updated to KiSAO 2.26 ([f25f243](https://github.com/biosimulations/Biosimulations/commit/f25f243bf85af086902badfd22b25547d974b9fb))
* **simulators:** added documentation for Python API conventions ([60669be](https://github.com/biosimulations/Biosimulations/commit/60669be5f166ac30a1af137240fbc049c13db331))

## [3.9.0](https://github.com/biosimulations/Biosimulations/compare/v3.8.0...v3.9.0) (2021-08-19)


### Features

* **dispatch:** added example run for MASSpy ([a589891](https://github.com/biosimulations/Biosimulations/commit/a589891a9367b58be7317f6b8f8da8545c2c44a7))
* **dispatch,ontology:** started to add MASS, RBA formats ([43a6153](https://github.com/biosimulations/Biosimulations/commit/43a615325f9695700ac2c9b68b2e124b0b03e3f9))
* **ontology:** updated to KiSAO 2.25 ([3fb5c54](https://github.com/biosimulations/Biosimulations/commit/3fb5c54af97d5631228c0e6456390377a867de5b))

## [3.8.0](https://github.com/biosimulations/Biosimulations/compare/v3.7.0...v3.8.0) (2021-08-18)


### Features

* **ontology:** updated to kisao 2.26 ([0f1f31a](https://github.com/biosimulations/Biosimulations/commit/0f1f31ae0383af38f9e7cd06aa28f022b7d6df07))

## [3.7.0](https://github.com/biosimulations/Biosimulations/compare/v3.6.0...v3.7.0) (2021-08-18)


### Bug Fixes

* update angular/cdk ([c7a18c1](https://github.com/biosimulations/Biosimulations/commit/c7a18c16f0c71fb0057a8cd331bbdd28abddb74b))
* **deps:** update dependency @openapi-contrib/openapi-schema-to-json-schema to v3.1.1 ([#2799](https://github.com/biosimulations/Biosimulations/issues/2799)) ([5d6453f](https://github.com/biosimulations/Biosimulations/commit/5d6453f1fc211b61ba7a0b40c9e18e7877bb874d))
* **deps:** update dependency @sendgrid/mail to v7.4.6 ([#2800](https://github.com/biosimulations/Biosimulations/issues/2800)) ([1e83398](https://github.com/biosimulations/Biosimulations/commit/1e83398a2b75f5af5b44506603f906f3c6e00b9d))
* **deps:** update dependency @typegoose/typegoose to v7.6.3 ([#2803](https://github.com/biosimulations/Biosimulations/issues/2803)) ([1a02d14](https://github.com/biosimulations/Biosimulations/commit/1a02d14610ec792e967f71cd4eec35c83af31d74))
* **deps:** update dependency auth0 to v2.36.1 ([#2823](https://github.com/biosimulations/Biosimulations/issues/2823)) ([fa11231](https://github.com/biosimulations/Biosimulations/commit/fa11231f161bb778e67f52287d707756b5315f20))
* **deps:** update dependency cache-manager to v3.4.4 ([#2804](https://github.com/biosimulations/Biosimulations/issues/2804)) ([3c62abb](https://github.com/biosimulations/Biosimulations/commit/3c62abb0576ede01ba5d9eec1ca0daf56a8b5a4e))
* **deps:** update nest ([#2807](https://github.com/biosimulations/Biosimulations/issues/2807)) ([875673f](https://github.com/biosimulations/Biosimulations/commit/875673f37fe2e4d0cc3f33497043f27899d03fd9))
* **dispatch-service:** added handling for case when no environment variables need to be set ([b51418d](https://github.com/biosimulations/Biosimulations/commit/b51418db31f38c534e774603070943d2abd0901e))
* **platform:** fix import of simulationrun metadata ([fb85650](https://github.com/biosimulations/Biosimulations/commit/fb85650554373a20e48be805d1068f0885bee7c0))


### Features

* **datamodel:** add common api query parameters ([c8bace5](https://github.com/biosimulations/Biosimulations/commit/c8bace5a78d6e991f0590a403927ab28757895f9))
* **dispatch:** added support for passing environment variables to simulators ([107221a](https://github.com/biosimulations/Biosimulations/commit/107221acfef0d611dba986f1da1e1782c1472d91))
* **dispatch-api:** add ability to get sparse simulationRuns ([5570ebb](https://github.com/biosimulations/Biosimulations/commit/5570ebb053df9d19c1fcb8f838564b1107d4c33f))
* **dispatch-api:** add endpoint for metadata ([24ffa40](https://github.com/biosimulations/Biosimulations/commit/24ffa40e96f42054ae4cee1f3eb62c6544ff08fe))
* **dispatch-api:** add tags, add ontology endpoint ([d056b4a](https://github.com/biosimulations/Biosimulations/commit/d056b4a74f0d19454e0eff0c515aea2eba5a0cd5))
* **dispatch-api:** create api model for metadata ([#2815](https://github.com/biosimulations/Biosimulations/issues/2815)) ([d55f7a5](https://github.com/biosimulations/Biosimulations/commit/d55f7a5d29e1f3cbb23ab472f84a2d8b961af843))
* **exceptions:** add better handling of validation errors ([e4e1986](https://github.com/biosimulations/Biosimulations/commit/e4e198623ae0d3f86460da18040faf9f28c9ee9a))

## [3.6.0](https://github.com/biosimulations/Biosimulations/compare/v3.5.0...v3.6.0) (2021-08-11)


### Bug Fixes

* **dispatch:** add logging when catching error ([85fab6c](https://github.com/biosimulations/Biosimulations/commit/85fab6cd89789cbf7357569194404e824bb60865))


### Features

* **dispatch:** added example runs for represillator model with SBML ([3400fa1](https://github.com/biosimulations/Biosimulations/commit/3400fa13eed1178ba74f76111b8ec83c995580f9))
* **dispatch:** added example simulation run for represillator model with OpenCOR ([ebffbae](https://github.com/biosimulations/Biosimulations/commit/ebffbae270595afd8d479d0a5ab6e90623b2323a))
* **dispatch:** added example simulation runs with visuaulizations using SBGN PD maps ([3c371e0](https://github.com/biosimulations/Biosimulations/commit/3c371e00758ab7a56b54afb5997da485fa3d071c))

## [3.5.0](https://github.com/biosimulations/Biosimulations/compare/v3.4.1...v3.5.0) (2021-08-09)


### Bug Fixes

* **deps:** pin dependency zone.js to 0.11.4 ([51bd9c7](https://github.com/biosimulations/Biosimulations/commit/51bd9c77245528ba8903d1644ddf5985899aa803))
* **dispatch:** added timeout for loading similar algorithms ([8c7a725](https://github.com/biosimulations/Biosimulations/commit/8c7a7256e13c9d6fd6c88fd4b5c89553d14906ef))
* **dispatch:** corrected display metadata loading indicator ([346aaeb](https://github.com/biosimulations/Biosimulations/commit/346aaeb809007d25b27e9a970f22eb7a3a715069))
* **dispatch:** fixed display of visualization loading indicator ([59fcd98](https://github.com/biosimulations/Biosimulations/commit/59fcd98686da9f8c6d22af7048b178d1b3362611))


### Features

* **dispatch:** added buttons to publish projects ([33c199a](https://github.com/biosimulations/Biosimulations/commit/33c199ae9e9ccdcc3ef41778f83c7caa7290d65e))
* **dispatch:** added display of errors with metadata of COMBINE archives ([41a0146](https://github.com/biosimulations/Biosimulations/commit/41a0146c8ac7421f13f0a9fe8aa61d2018f38ff0))
* **dispatch:** added support for XPP ([3d108ad](https://github.com/biosimulations/Biosimulations/commit/3d108adb5e468d408b4fbad061eab48f84861ea0))
* **dispatch:** improving capabilities when COMBINE service is down ([dc9ca09](https://github.com/biosimulations/Biosimulations/commit/dc9ca095ff96405e46a13d5f7ce56b95645308e7))
* **dispatch:** simplified designing 2D line/scatter plots with multiple curves ([beab68e](https://github.com/biosimulations/Biosimulations/commit/beab68e194e637940c57637520f292d0262fa328))

## [3.4.1](https://github.com/biosimulations/Biosimulations/compare/v3.4.0...v3.4.1) (2021-07-30)


### Bug Fixes

* **deps:** pin dependency @ngbmodule/material-carousel to 0.7.1 ([0d1acde](https://github.com/biosimulations/Biosimulations/commit/0d1acde594b0bc455c50209603684b8da7f66a02))
* **dispatch-api:** send correct message when simulation status changes ([a3c9c62](https://github.com/biosimulations/Biosimulations/commit/a3c9c6235102a33dafbc8414d0d5535c1a641f2f)), closes [#2739](https://github.com/biosimulations/Biosimulations/issues/2739)

## [3.4.0](https://github.com/biosimulations/Biosimulations/compare/v3.3.0...v3.4.0) (2021-07-29)


### Bug Fixes

* **dispatch:** corrected processing of metadata while status is pinging ([e29fb0e](https://github.com/biosimulations/Biosimulations/commit/e29fb0e54343da3c43db5f949dc48e83842845ca))


### Features

* **dispatch:** added example simulation run for activity flow diagram ([f534c33](https://github.com/biosimulations/Biosimulations/commit/f534c33835274eb113c9f62209d13810cb55f778))
* **dispatch:** expanded support for connecting SED-ML to Vega ([439bbeb](https://github.com/biosimulations/Biosimulations/commit/439bbebcfcc562304018711b9e0e485cba099eb1))
* **ontology:** updated SBO for additional framework terms ([5eaa097](https://github.com/biosimulations/Biosimulations/commit/5eaa097753353e9134a81b92554f3ca7efd8335e))

## [3.3.0](https://github.com/biosimulations/Biosimulations/compare/v3.2.0...v3.3.0) (2021-07-23)


### Bug Fixes

* **dispatch:** hiding figures/tables section when there are no figures/tables ([056caf8](https://github.com/biosimulations/Biosimulations/commit/056caf8d6f7ec6c7b0f6ad2116588e8f6dea751d))


### Features

* **dispatch,simulators:** added documentation about generating data visualizations ([0066522](https://github.com/biosimulations/Biosimulations/commit/00665225875e75a62a93dd93fd545f8e823f9ecc))
* making creation data metadata optional ([7812d65](https://github.com/biosimulations/Biosimulations/commit/7812d654f1ddfed9f9c2ea00b63ae31c3a537942))

## [3.2.0](https://github.com/biosimulations/Biosimulations/compare/v3.1.0...v3.2.0) (2021-07-22)


### Bug Fixes

* **dispatch-api:** change path from 'run' to 'runs' ([ead8d80](https://github.com/biosimulations/Biosimulations/commit/ead8d807ffe48a2cb50e54d88a65e880d00b6a70))


### Features

* **dispatch:** improved Vega error handling ([56a1e0c](https://github.com/biosimulations/Biosimulations/commit/56a1e0ce1f7147130bd18651dcac9d0b6953bb09))
* **dispatch:** updated example runs for new vis and metadata ([1683451](https://github.com/biosimulations/Biosimulations/commit/168345199fb240e098f168211609973076251a0b))
* **platform,platform-api:** platform gets projects from api ([f0b010d](https://github.com/biosimulations/Biosimulations/commit/f0b010d68b592765acb172c27a1b527ca4d9d157))
* **simulators:** added documentation about recommendation to use Identifiers.org URIs ([5445b08](https://github.com/biosimulations/Biosimulations/commit/5445b08a46678b4f73ea25fec53b38c9fdc6de4d))

## [3.1.0](https://github.com/biosimulations/Biosimulations/compare/v3.0.2...v3.1.0) (2021-07-19)


### Bug Fixes

* **deps:** pin dependencies ([5c762b6](https://github.com/biosimulations/Biosimulations/commit/5c762b64117e22863f8edf96c0d358256536b765))
* **deps:** update nest monorepo ([790aa52](https://github.com/biosimulations/Biosimulations/commit/790aa52b225b7d1eca6312da65fa0ff6b3d6fb9c))
* **dispatch:** fixed handling of query arguments to run simulation route ([546c49b](https://github.com/biosimulations/Biosimulations/commit/546c49bc344dfade3fd677998a665066e35ea2ce))
* **simulators:** reenabling display of simulator validation test results ([902bc93](https://github.com/biosimulations/Biosimulations/commit/902bc9374875572f4072ffde72d369a60f6f31c1)), closes [#2696](https://github.com/biosimulations/Biosimulations/issues/2696)


### Features

* **dispatch:** added 1D histogram plot ([379e0a7](https://github.com/biosimulations/Biosimulations/commit/379e0a7d74ad74b1cdc35e7babba9be7b1155c01))
* **dispatch:** added 2D heatmap data visualization ([147b6ad](https://github.com/biosimulations/Biosimulations/commit/147b6ad4ae52e77b4c02eac6cfc5b17dacf8a5c1))
* **dispatch:** added ability to add files to COMBINE archives ([c05631e](https://github.com/biosimulations/Biosimulations/commit/c05631e52fbcba6c72ccb24380b2d66167790218))
* **dispatch:** added exporting user-configured histogram viz to Vega ([6342e44](https://github.com/biosimulations/Biosimulations/commit/6342e448059fcc100d37b7246206cfb01b1bd8e8))
* **dispatch:** added Vega export for 2D heatmap, improved visualization form validation ([784917f](https://github.com/biosimulations/Biosimulations/commit/784917f3172e2416fe3e959218ec27e452fe4a79))
* **dispatch:** added Vega export for 2d line/scatter plot ([8f2cff0](https://github.com/biosimulations/Biosimulations/commit/8f2cff046bc3db7fdab6a8690a780da1d3ae7867))
* **dispatch:** improved plotting ([852545b](https://github.com/biosimulations/Biosimulations/commit/852545b75fb4dd1bca2b0594ca519f1cab9a111d))
* **dispatch:** linked Vega signals to attributes of SED-ML simulations ([88a68c3](https://github.com/biosimulations/Biosimulations/commit/88a68c332c9f630bc9c5cffa02b9c9e2e3f0058a))
* **dispatch:** linking published figures/tables to displayed visualizations ([44c810a](https://github.com/biosimulations/Biosimulations/commit/44c810acdc8ed31173387bd5521dbc03c093008a))
* **platform:** implement viewing a project ([0ad9af3](https://github.com/biosimulations/Biosimulations/commit/0ad9af33635272cfc29e4cd9b3d2d71cdb03dbe4))
* **platform-api:** add skeleton implementation ([0758052](https://github.com/biosimulations/Biosimulations/commit/0758052e88b185d72e3cedc46256377a8e3d9753))

## [3.0.2](https://github.com/biosimulations/Biosimulations/compare/v3.0.1...v3.0.2) (2021-07-13)


### Bug Fixes

* **mail-service,dispatch-service:** fix nats-server connection ([9143d0c](https://github.com/biosimulations/Biosimulations/commit/9143d0cabba6a22d44398cd4eec1cb1a5033e2f9))

## [3.0.1](https://github.com/biosimulations/Biosimulations/compare/v3.0.0...v3.0.1) (2021-07-13)


### Bug Fixes

* try new server options ([1851742](https://github.com/biosimulations/Biosimulations/commit/1851742fdb222a6330a6b2f52b814ee7d3273c5a))

## [3.0.0](https://github.com/biosimulations/Biosimulations/compare/v2.5.2...v3.0.0) (2021-07-13)


### Bug Fixes

* **mail-service,dispatch-service:** fix import of http module ([805e48f](https://github.com/biosimulations/Biosimulations/commit/805e48f34dec1b50308c20affa786fac1ae646f5))


### Features

* **simulators-api:** add a query argument to include the results of the validation ([710be08](https://github.com/biosimulations/Biosimulations/commit/710be085ec732d851aa89e78773c4ba12e7e682e)), closes [#2668](https://github.com/biosimulations/Biosimulations/issues/2668)


### BREAKING CHANGES

* **simulators-api:** validation data is no longer returned by default. A Query argument is needed to include the validation information

## [2.5.2](https://github.com/biosimulations/Biosimulations/compare/v2.5.1...v2.5.2) (2021-07-13)


### Bug Fixes

* type and build fixes ([6812bd0](https://github.com/biosimulations/Biosimulations/commit/6812bd0de12c3716ac8928154ed73aa92953dc40))
* **dispatch-api:** fix error with parsing outputIds ([9fac99f](https://github.com/biosimulations/Biosimulations/commit/9fac99f68094919de23d632795bf131b9fb8a1ef)), closes [#2683](https://github.com/biosimulations/Biosimulations/issues/2683)

## [2.5.1](https://github.com/biosimulations/Biosimulations/compare/v2.5.0...v2.5.1) (2021-07-09)


### Bug Fixes

* **simulators-api:** Allow for date based versions ([0c8fb8d](https://github.com/biosimulations/Biosimulations/commit/0c8fb8d00a36675f652a665d9c279709e798c212)), closes [#2681](https://github.com/biosimulations/Biosimulations/issues/2681)

## [2.5.0](https://github.com/biosimulations/Biosimulations/compare/v2.4.0...v2.5.0) (2021-07-09)


### Bug Fixes

* **dispatch:** downloading created COMBINE/OMEX archives ([03895bf](https://github.com/biosimulations/Biosimulations/commit/03895bf82abfa0b6d9c7c7db186a31d29e726b49))
* correcting size of form fields ([ae69630](https://github.com/biosimulations/Biosimulations/commit/ae69630cf3e45abfbcaaec17787956e7189e5e53))
* **shared-exceptions:** include error metadata in the "meta" output ([2be0178](https://github.com/biosimulations/Biosimulations/commit/2be0178af9003156ad25fa22e8c7fe51457c9556))


### Features

* **combine-service:** adding support for creating steady-state analyses of logical models ([a9e6667](https://github.com/biosimulations/Biosimulations/commit/a9e6667034c2b35d4379dff72a2d7cefe4d4f4d8))
* **combine-service:** updating to biosimulators-utils 0.1.93 ([ca0a21e](https://github.com/biosimulations/Biosimulations/commit/ca0a21e33d7c2a54f8bd6d9aa9d8c6943da955b2))
* **shared-exceptions:** add validation pipe error factory ([35edb4d](https://github.com/biosimulations/Biosimulations/commit/35edb4d4f73d82e4bbb17bbd701d13fc580093af))

## [2.4.0](https://github.com/biosimulations/Biosimulations/compare/v2.3.0...v2.4.0) (2021-07-08)


### Bug Fixes

* **combine-service:** add protocol to server in API spec ([bcf4119](https://github.com/biosimulations/Biosimulations/commit/bcf41192f15894f7993239f16b14f415c1a85910))
* **combine-service:** debugged specifications ([f8b9420](https://github.com/biosimulations/Biosimulations/commit/f8b9420fadbcf10c367f26c6b77bf736dc5463f3))
* **combine-service:** fix api spec ([ebab457](https://github.com/biosimulations/Biosimulations/commit/ebab457a262b594744e81cfdfe4bfdd9af45a4c2))
* **platform:** enable strict template checking, fix type errors ([9facdb1](https://github.com/biosimulations/Biosimulations/commit/9facdb19c1cacae3d22d4bd952c0f1f0cbaf8035)), closes [#2185](https://github.com/biosimulations/Biosimulations/issues/2185)


### Features

* **combine-service:** add combine-service api client library ([bfa25b8](https://github.com/biosimulations/Biosimulations/commit/bfa25b8f0319de62d3c6a5902597d51c68b8eb96))
* **dispatch:** added example simulation run for GINsim ([3e639b3](https://github.com/biosimulations/Biosimulations/commit/3e639b30c49e3149d966518d5a82908a3905e831))
* **dispatch:** added example simulation runs for GINsim, LibSBMLSim ([06ad0a4](https://github.com/biosimulations/Biosimulations/commit/06ad0a4be78dc595305492e975ba856722d27b27))
* **dispatch:** adding support for GINML, ZGINML to COMBINE archive creation and execution ([9f949e5](https://github.com/biosimulations/Biosimulations/commit/9f949e561256288f2851468a83c96160cc14f7fe))
* **dispatch,ontology:** add terms for GINsim format ([22d8a7b](https://github.com/biosimulations/Biosimulations/commit/22d8a7b2daec93c086bacc3c61a279dc85481cfd))
* **ontology:** updating to KiSAO 2.19 with terms for logical modeling ([c47e63b](https://github.com/biosimulations/Biosimulations/commit/c47e63b7767eecde90c033d8c11ef55b89678d4a))
* **ontology,combine-service:** update to KiSAO 2.20 ([fadf3da](https://github.com/biosimulations/Biosimulations/commit/fadf3da7c0e714267a89ac903acf516be5f00533))

## [2.3.0](https://github.com/biosimulations/Biosimulations/compare/v2.2.1...v2.3.0) (2021-07-02)

### Features

- **dispatch:** Added tab to simulation run page to display metadata about the simulation project ([#2667](https://github.com/biosimulations/Biosimulations/issues/2667)) ([dde87fa](https://github.com/biosimulations/Biosimulations/commit/dde87faae5e558c3bbe86f6f17467ae747da55d8)), closes [#2661](https://github.com/biosimulations/Biosimulations/issues/2661)

## [2.2.1](https://github.com/biosimulations/Biosimulations/compare/v2.2.0...v2.2.1) (2021-07-01)

### Bug Fixes

- **dispatch:** fix example simulation runs ([60d91c1](https://github.com/biosimulations/Biosimulations/commit/60d91c1bb70e6ae08274a9380143baa19fa51043)), closes [#2653](https://github.com/biosimulations/Biosimulations/issues/2653)
- **simulators-api:** fix getting latest version ([4594c96](https://github.com/biosimulations/Biosimulations/commit/4594c96b53859e03960458cd001cf8614d64f64c)), closes [#2664](https://github.com/biosimulations/Biosimulations/issues/2664)

## [2.2.0](https://github.com/biosimulations/Biosimulations/compare/v2.1.0...v2.2.0) (2021-06-30)

### Bug Fixes

- **dispatch:** correct integration between simulation results and SED plots ([0bab60f](https://github.com/biosimulations/Biosimulations/commit/0bab60fe06cc52d55a670d8957e385dc7f247854))
- **dispatch:** download file instead of redirect ([cd2840d](https://github.com/biosimulations/Biosimulations/commit/cd2840d98d84f13eab34cea479a09da23187fe14)), closes [#2435](https://github.com/biosimulations/Biosimulations/issues/2435)
- **dispatch:** use correct api to get simulator info ([1e66f1f](https://github.com/biosimulations/Biosimulations/commit/1e66f1f85f4436987ca034c3cdafad9536c12b9e))

### Features

- add some shared endpoints ([567e4c2](https://github.com/biosimulations/Biosimulations/commit/567e4c27de05655d3b78b441e84231977afd234b))
- **auth-client:** Cache tokens locally ([f53c9f8](https://github.com/biosimulations/Biosimulations/commit/f53c9f8d4c9c3e2bed497ec85c4c53d774af9fb1)), closes [#2503](https://github.com/biosimulations/Biosimulations/issues/2503)
- **auth-common:** add util functions ([e0ac842](https://github.com/biosimulations/Biosimulations/commit/e0ac842518af8e6909493cb1b2b774a56faf6b17))

### Performance Improvements

- **dispatch-service:** use /local as the singularity cache/working directory ([c63b58c](https://github.com/biosimulations/Biosimulations/commit/c63b58c35a0c3da71910523a3baf0f445f5e493a))

## [2.1.0](https://github.com/biosimulations/Biosimulations/compare/v2.0.0...v2.1.0) (2021-06-18)

### Bug Fixes

- **dispatch-service:** remove check for process flag ([f7f88cc](https://github.com/biosimulations/Biosimulations/commit/f7f88cce2fbc54df13e34ef5212f1491036ec8b5)), closes [#2577](https://github.com/biosimulations/Biosimulations/issues/2577)

### Features

- **dispatch-api, dispatch-service:** add status reason to datamodel ([ca9bcb6](https://github.com/biosimulations/Biosimulations/commit/ca9bcb6c7d7ffcb0328ef679d5a82801995add45)), closes [#2441](https://github.com/biosimulations/Biosimulations/issues/2441)

## [2.0.0](https://github.com/biosimulations/Biosimulations/compare/v1.0.0...v2.0.0) (2021-06-17)

### Bug Fixes

- **dispatch-api:** bind class to this variable in map ([b4bb3ca](https://github.com/biosimulations/Biosimulations/commit/b4bb3ca27cd52d27abe68dcaa524a158a1a73507))
- dispatch frontend uses the updated api parameter ([#2636](https://github.com/biosimulations/Biosimulations/issues/2636)) ([a13779c](https://github.com/biosimulations/Biosimulations/commit/a13779cdc320d58c595f85399ca4d7747d603657)), closes [#2635](https://github.com/biosimulations/Biosimulations/issues/2635)

### Features

- **dispatch-api, dispatch-service:** Use HSDS to get simulation run data ([33b8030](https://github.com/biosimulations/Biosimulations/commit/33b8030e60fcbd2eb693e2a962620cf42855b4e4)), closes [#2533](https://github.com/biosimulations/Biosimulations/issues/2533) [#2442](https://github.com/biosimulations/Biosimulations/issues/2442) [#2440](https://github.com/biosimulations/Biosimulations/issues/2440) [#2369](https://github.com/biosimulations/Biosimulations/issues/2369) [#2069](https://github.com/biosimulations/Biosimulations/issues/2069)

### BREAKING CHANGES

- **dispatch-api, dispatch-service:** Dispatch API no longer has endpoints for creating or updating "Result" objects.
  The output of the results endpoints are updated to include information about type and shape of the data.
  The parameter "sparse" has been changed to "includeData".
  The datamodel for results has been adjusted to include all outputs, not just reports. "reports" has been renamed to "outputs"

# 1.0.0 (2021-06-16)

This is an arbitrary starting point for tracking changes and versioning. It should not be considered as the "first release".

### Bug Fixes

- bash script ([866b58a](https://github.com/biosimulations/Biosimulations/commit/866b58a244d3483fa6afa2ae0e8383e234920ba6))
- add check for large files downloading ([ca10aa5](https://github.com/biosimulations/Biosimulations/commit/ca10aa5f2c44fafcff8471c0da0bbd9db2a655ed)), closes [#2536](https://github.com/biosimulations/Biosimulations/issues/2536)
- bring inline with datamodel ([af53a54](https://github.com/biosimulations/Biosimulations/commit/af53a54a5e5bc91835e114e67e10fc69883c7f9b))
- change url to download results ([3d264d4](https://github.com/biosimulations/Biosimulations/commit/3d264d4abfdb85aa4cae99e1477cbf4666b8ba36)), closes [#2561](https://github.com/biosimulations/Biosimulations/issues/2561)
- check job status after completion ([c16649c](https://github.com/biosimulations/Biosimulations/commit/c16649c507c35f1e086d41fb496a573549f925ba))
- cleanup logs ([fbd330f](https://github.com/biosimulations/Biosimulations/commit/fbd330f2762f6192b7e269ffcd2e8ede93c5ad14))
- correct value for constant ([e2d3a68](https://github.com/biosimulations/Biosimulations/commit/e2d3a68a39f5e6ad3216daf109087a2b1c43f26b))
- fix error in reading port ([e1f6fb9](https://github.com/biosimulations/Biosimulations/commit/e1f6fb923a42283d1b42765b4d0376a146f406ef))
- fix logs and context buttons ([777e8e8](https://github.com/biosimulations/Biosimulations/commit/777e8e8f79f829b3762c7aa189a9d6184f4b24a1)), closes [#2543](https://github.com/biosimulations/Biosimulations/issues/2543) [#2540](https://github.com/biosimulations/Biosimulations/issues/2540)
- fix redis queue and port ([5f33a19](https://github.com/biosimulations/Biosimulations/commit/5f33a192203323e30d6badd4b6500cc056b3ef34))
- fix s3 key for downloading outputs ([f585a9a](https://github.com/biosimulations/Biosimulations/commit/f585a9a6295cbbf60e73c05d5ae908713d1ef5ee)), closes [#2622](https://github.com/biosimulations/Biosimulations/issues/2622)
- fix spelling of library ([a471e95](https://github.com/biosimulations/Biosimulations/commit/a471e95e6684ee093d036f343efbba50df327563))
- fix test ([6f236df](https://github.com/biosimulations/Biosimulations/commit/6f236df6b5186e44ccf9c459faa83698cf22d7ae))
- fix test ([6af0ca8](https://github.com/biosimulations/Biosimulations/commit/6af0ca8a0b9f2d557a0fd416475151261a46fb88))
- lint fix ([a26c24b](https://github.com/biosimulations/Biosimulations/commit/a26c24b17e9d1f72bdb53860b7f27b300030ec68))
- order of operations for creating results ([eac31e0](https://github.com/biosimulations/Biosimulations/commit/eac31e01bde327b1d3ff89f6a7cf7480e5d0c96d))
- propely set name and filetype of outputs ([951c239](https://github.com/biosimulations/Biosimulations/commit/951c239983ea93009565287a8ac9bfd3deae8052))
- Remove bad library import ([ecc86fa](https://github.com/biosimulations/Biosimulations/commit/ecc86fa6d9abf59b466ea02d25d62c1119d07de8)), closes [#2420](https://github.com/biosimulations/Biosimulations/issues/2420)
- remove xdg runtime directory ([f5ec15b](https://github.com/biosimulations/Biosimulations/commit/f5ec15bd726ab4afa01b0c2be4688217d4d89198))
- resolve build errors ([6691ebe](https://github.com/biosimulations/Biosimulations/commit/6691ebedbda107862cbf731cb891044c426e5fc9))
- typo in return statement ([1f6c4fc](https://github.com/biosimulations/Biosimulations/commit/1f6c4fc0bda0780170530789e27e4fab4233d2e3))
- update default stoage URL ([f9b0d75](https://github.com/biosimulations/Biosimulations/commit/f9b0d75b3c4370653d1c0596157cc381fa6573f0))
- update logs ([818a0c3](https://github.com/biosimulations/Biosimulations/commit/818a0c347529c42d697ac972c15e17c09b5e0372))
- update sbatch memoy amount ([b9026f9](https://github.com/biosimulations/Biosimulations/commit/b9026f96ff2e5b4876559b1b116c1d5cdebbfb8d))
- update sbatch script to use custom module ([0ef1c52](https://github.com/biosimulations/Biosimulations/commit/0ef1c52de4d6703032decaff9b2c8941175c70fb))
- use job status to determine completion ([adb12a0](https://github.com/biosimulations/Biosimulations/commit/adb12a0efbe07e82346ddada18ad93342d1cede5))
- **apps/frontend:** relative import ([3854f27](https://github.com/biosimulations/Biosimulations/commit/3854f272fdd21847e522cb03f25353b06a3c3028))
- **auth:** check for logged in before intercepting ([7c22a19](https://github.com/biosimulations/Biosimulations/commit/7c22a19a2a33cff63067d30dd19ff8bfe091a189))
- **auth:** Check for username correctly ([05996b3](https://github.com/biosimulations/Biosimulations/commit/05996b376ed7b08a4de974cf39029cb9956c1070))
- **disatch:** Fix open api schema ([ff15503](https://github.com/biosimulations/Biosimulations/commit/ff15503862e62412418f3144353e76a8cd877f7b))
- **dispatch:** Fix observable piping ([34e7086](https://github.com/biosimulations/Biosimulations/commit/34e7086261b20ffdf42b6edcf3112f972432ea79))
- **dispatch:** parsing results accounts for quotes ([d055a7b](https://github.com/biosimulations/Biosimulations/commit/d055a7ba328e1e4d4d5f8094a661377a8e5294f9)), closes [#2459](https://github.com/biosimulations/Biosimulations/issues/2459)
- **dispatch:** patch error handling ([d2d98e5](https://github.com/biosimulations/Biosimulations/commit/d2d98e57bd1632d8f289e2d9fd017443a653c8db))
- **dispatch:** remove bad environment variables ([3c31b7d](https://github.com/biosimulations/Biosimulations/commit/3c31b7de39d05b5bfaf640a1553a204267247eab)), closes [#2476](https://github.com/biosimulations/Biosimulations/issues/2476)
- **dispatch:** Simulation results not saved for some simulations and overall status doesn't reflect such errors ([#2428](https://github.com/biosimulations/Biosimulations/issues/2428)) ([acd2dff](https://github.com/biosimulations/Biosimulations/commit/acd2dff837834e6732f4b5074c433f90a9523d06)), closes [#2416](https://github.com/biosimulations/Biosimulations/issues/2416)
- use https for auth0 image ([19a4dcc](https://github.com/biosimulations/Biosimulations/commit/19a4dcc53d6572be6b60a9bb8a4d9db4bd89afc6))
- **forms:** connect taxon form properly ([a4088f7](https://github.com/biosimulations/Biosimulations/commit/a4088f79bf3a6c06cfc0dfbb48936d820428c378))
- **forms:** fix reference form implementation ([5b3eba4](https://github.com/biosimulations/Biosimulations/commit/5b3eba498a464fe77f7bbad3ae3365db8bb3e6bc))
- **forms:** fix some taxon form details ([8802ced](https://github.com/biosimulations/Biosimulations/commit/8802cedb84fc7c1afdc2b5cea244bc11c5a67399))
- **forms:** fix tags form ([b60b99c](https://github.com/biosimulations/Biosimulations/commit/b60b99cc6c8526a7403c585f8f79315c757972a5))
- **forms:** fix validation error with taxon form ([6cb2d43](https://github.com/biosimulations/Biosimulations/commit/6cb2d43e7edda5edff8310cc5bcffc5a95cedec8))
- **forms:** Form does not scroll over topbar ([70f3275](https://github.com/biosimulations/Biosimulations/commit/70f327598c2cacd278625dd78c16e93761974d54))
- **forms:** set file form disable properly ([c523d0e](https://github.com/biosimulations/Biosimulations/commit/c523d0e8230c14c3c64a0266107b7777851ae114))
- **gaurds:** Gaurd loads underconstruction pages by default ([9f7a810](https://github.com/biosimulations/Biosimulations/commit/9f7a810bd6aa000039fb7dd3e3e5b421f63118e3))
- **grid:** grid work with async resources ([d621cd3](https://github.com/biosimulations/Biosimulations/commit/d621cd3d11402aca909a005e376d214160f117b2))
- **interceptor:** API token, error handling ([78be137](https://github.com/biosimulations/Biosimulations/commit/78be1377afbf3877885744074fe8dccd0ffe1ca6))
- **interceptor:** Fixed a bug in the error handling of the interceptor ([1397c96](https://github.com/biosimulations/Biosimulations/commit/1397c966e4a83fd1e0b660daeef092617cecc106))
- **models:** Edit component calls subscribe ([50f15ac](https://github.com/biosimulations/Biosimulations/commit/50f15accb120fcdb136262541e386bacee630cb3))
- **navigation:** get username via async ([7415eba](https://github.com/biosimulations/Biosimulations/commit/7415ebae27841b6147931820869b8ef06a55cc89))
- **navigation:** have navigation work with async ([0da33ea](https://github.com/biosimulations/Biosimulations/commit/0da33eaf426d907115ba07636674261611e73a7c))
- **polyfill:** add back polyfill ([30e9c1d](https://github.com/biosimulations/Biosimulations/commit/30e9c1de6204f3bba371be0ad53c89b8a8939f1f))
- **resource service:** add query params to read ([d1d7d38](https://github.com/biosimulations/Biosimulations/commit/d1d7d38610582ed45769b5f3d404f9222cbd08cc))
- **resources:** Move more functionality into abstract class ([5e94323](https://github.com/biosimulations/Biosimulations/commit/5e94323cdf7405de6d80b69bb18a7217d02e9922))
- **serializer:** fix private public being flipped ([035241f](https://github.com/biosimulations/Biosimulations/commit/035241f949d25c38bf9767a85ee805df958d0395))
- **serializers:** improve serializers ([0d3f004](https://github.com/biosimulations/Biosimulations/commit/0d3f004175ece11a48c668959d9b5b3e217994ae))
- **serializers:** user serializer returns none for '' ([221ce41](https://github.com/biosimulations/Biosimulations/commit/221ce413570b04671da1e77581d441a423eb4647))
- **services:** Dont cast http reponse to resource ([b20e69d](https://github.com/biosimulations/Biosimulations/commit/b20e69d45dcdcaf69dedc41cc157dc9cfee14005))
- **simulations:** fix async view ([d2cad0f](https://github.com/biosimulations/Biosimulations/commit/d2cad0fd3230714bae9be8858389456bcb6ac9ca))
- **tests:** fix common test issues ([61c7621](https://github.com/biosimulations/Biosimulations/commit/61c76219e91b7b43675014ea7293dd20e3dc1cf2))
- **tests:** fix common test issues ([0c4e5f3](https://github.com/biosimulations/Biosimulations/commit/0c4e5f3f9fa58cc600a45b401a5a4ab3ee23c315))
- **user:** profile edit component properly creates user ([0bed682](https://github.com/biosimulations/Biosimulations/commit/0bed68288d6fcb1c885138ec83e0c7309f3415d3))
- **visualizations:** fix licence view ([f5d9973](https://github.com/biosimulations/Biosimulations/commit/f5d997394786baf2fac0fe4a0b6ea11c402b9881))

### Features

- add client library ([820647b](https://github.com/biosimulations/Biosimulations/commit/820647b4c13ddae30174232c1da0cd5f88990dc6))
- add config for queue ([c7ec4a1](https://github.com/biosimulations/Biosimulations/commit/c7ec4a1ed6b840e11968597f5c530bf6d0a15566))
- add hsds client module ([1730521](https://github.com/biosimulations/Biosimulations/commit/17305212950d8f27e3d577660de00828b90a5f2d))
- use org for getting latest simulator ([2f7c503](https://github.com/biosimulations/Biosimulations/commit/2f7c503b30f921f8e21919934774c13dec4113d8))
- **api:** add config to api ([931dcf5](https://github.com/biosimulations/Biosimulations/commit/931dcf518c7d57adaeb767bdc18e7fcee151be1a))
- **api:** add crud skeleton for routes ([64fce18](https://github.com/biosimulations/Biosimulations/commit/64fce18062fa43ea5c3063fe29a569d8b97f1d09))
- **api:** add open api spec generation ([659f8b4](https://github.com/biosimulations/Biosimulations/commit/659f8b4a51cfb201aea00fc7d5b2bfac6c7a9c14))
- **author form:** build out author form ([afe666a](https://github.com/biosimulations/Biosimulations/commit/afe666ae014b6463d53349913488d894f64f5aa6))
- **datamodel:** add gaurds ([2ffdf04](https://github.com/biosimulations/Biosimulations/commit/2ffdf04605f2ecb9dd41d7c7c681e5843ece23c2))
- **datamodel:** add more of the core datamodel ([9b28f83](https://github.com/biosimulations/Biosimulations/commit/9b28f834cc0a671e1628617d97acc6bae9068412))
- **datamodel:** add properties to format ([67b024c](https://github.com/biosimulations/Biosimulations/commit/67b024ccb1620937a711da7486bc328f9f698328)), closes [#462](https://github.com/biosimulations/Biosimulations/issues/462)
- **datamodel:** add url to ontology ([3c8a169](https://github.com/biosimulations/Biosimulations/commit/3c8a1699b279acd687b1fcd80a420bb509cd09fe))
- **datamodel:** redefine core objects as set of attriutes ([95ea0a7](https://github.com/biosimulations/Biosimulations/commit/95ea0a78c2c23617b8ed103de55e41fdf2059122))
- **datamodel:** redfine core resources as primary and secondary ([dd2fec4](https://github.com/biosimulations/Biosimulations/commit/dd2fec4ac6d0f8fce9bfa3d1232e95b1dd528986))
- **dispatch:** add a dispatch service ([152b3b0](https://github.com/biosimulations/Biosimulations/commit/152b3b06536428850af601187eaf7f243f45b4d6))
- **errors:** add default errors component ([a10b927](https://github.com/biosimulations/Biosimulations/commit/a10b92773d263c289e59bab261381fb49bf3f953))
- **errors:** Create 404 component ([e1c7d33](https://github.com/biosimulations/Biosimulations/commit/e1c7d334ccde6684284ba1cab0ce38b2f109c3e9))
- **errors:** Create errors module ([7e3cac4](https://github.com/biosimulations/Biosimulations/commit/7e3cac4e27b4104595f94ef269e478916336b168))
- **errors:** Create under construction component ([5d2fbe1](https://github.com/biosimulations/Biosimulations/commit/5d2fbe1405117f8acc52b7ae5cbc9c5b7d4f5974))
- **errors:** slight changes to underConstruction ([bb06012](https://github.com/biosimulations/Biosimulations/commit/bb06012172c2a86c014ad1176832fe68a3439428))
- **forms:** add a component for file inputs ([e843f87](https://github.com/biosimulations/Biosimulations/commit/e843f873114fde51c917d01f48f6ce3f3fe306ff))
- **forms:** add abstract array subform ([1d70e47](https://github.com/biosimulations/Biosimulations/commit/1d70e479618af355416e287f973b5b8bbdc3db67))
- **forms:** add authors and identifiers components ([d2e35e8](https://github.com/biosimulations/Biosimulations/commit/d2e35e8bca295aa5a536eea3fbabc02c9119c13a))
- **forms:** Add identifier form control ([dd80d37](https://github.com/biosimulations/Biosimulations/commit/dd80d37d640c5eae612d216f26e3472ab5607cd8))
- **forms:** add name form control ([8557559](https://github.com/biosimulations/Biosimulations/commit/8557559e2e2d85e0409e31aa9cb2faaff39165a3))
- **forms:** add required validator to fields ([ceaacb5](https://github.com/biosimulations/Biosimulations/commit/ceaacb53174505e27ab2e98da1d80db4d346a6bf))
- **forms:** add resource form skeleton ([e7a56d5](https://github.com/biosimulations/Biosimulations/commit/e7a56d5892c80205fba48809e1f89dc347f3abff))
- **forms:** add taxon form ([48a2e86](https://github.com/biosimulations/Biosimulations/commit/48a2e8605a3de64f75368fd8a4bd0bcd8b6d6785))
- **forms:** add to resource form implementation ([dc0b9ed](https://github.com/biosimulations/Biosimulations/commit/dc0b9ed74e1ec8435e6422d8623c0b5aa75344b4))
- **forms:** add username form ([7ba4e81](https://github.com/biosimulations/Biosimulations/commit/7ba4e814ca1b5155bd5ae067c9385dbe3e84a30c))
- **forms:** create a model form ([353d2ab](https://github.com/biosimulations/Biosimulations/commit/353d2ab55ca5c4d6d7a3539a55e09604888bc5ef))
- **forms:** Create descriptions form control ([87b0fdd](https://github.com/biosimulations/Biosimulations/commit/87b0fdd95967957dbcc2033e858df9f975d838cb))
- **forms:** Create edit-preview component ([5177f05](https://github.com/biosimulations/Biosimulations/commit/5177f05ee5d86c6f19ee735b9c181f939971150d))
- **forms:** create model format form ([af0ad55](https://github.com/biosimulations/Biosimulations/commit/af0ad55b1972146aa6313e89dd8359a7d5d6d604))
- **forms:** enable access form component ([741afa3](https://github.com/biosimulations/Biosimulations/commit/741afa323769ffc6209e1b7cd6056311e3625a65))
- **forms:** Enable Drag/Drop ([799147b](https://github.com/biosimulations/Biosimulations/commit/799147b2277a7cc06191a50eb3cd6f66cc589f13))
- **forms:** Finalize author form ([3f4f86c](https://github.com/biosimulations/Biosimulations/commit/3f4f86cfb56c7d3785362f384a1a3378310a9cfd))
- **forms:** Generalize single field controls ([0da8122](https://github.com/biosimulations/Biosimulations/commit/0da8122c33d69a36a0a8db810e60a40a9340404d))
- **forms:** implement licence form ([2d2d9b1](https://github.com/biosimulations/Biosimulations/commit/2d2d9b1340e30867321ce8fbcef86c1e7e3be8e4))
- **forms:** implement refrences form ([122fa33](https://github.com/biosimulations/Biosimulations/commit/122fa33703cbafcaef32a4665ad130ee8599cd29))
- **forms:** implement resource form ([85e89ea](https://github.com/biosimulations/Biosimulations/commit/85e89ea97e8990e23b9b9477e9017b3f289854a1))
- **forms:** implement tags form ([0b1e480](https://github.com/biosimulations/Biosimulations/commit/0b1e48018e8cc1f22620e775cbd836fb63e315bb))
- **forms:** improve disable handling ([88cace8](https://github.com/biosimulations/Biosimulations/commit/88cace866e81e97def3dbbbedeacaf4daf24f435))
- **forms:** styling ([32044ef](https://github.com/biosimulations/Biosimulations/commit/32044ef922ae05010aebe1c9e83016007e7b069d))
- **home:** Add sponsors section ([8f84a67](https://github.com/biosimulations/Biosimulations/commit/8f84a677f550655b17d425dc6c961276a9fcaca5))
- **logging:** Added logging ([18476ca](https://github.com/biosimulations/Biosimulations/commit/18476ca216a6f0ddd5274d51073e7733c2be0c15))
- **login:** add styling ([642a71c](https://github.com/biosimulations/Biosimulations/commit/642a71cb269498b8272dd0a0b0e85c2f395fa170))
- **login:** redirect works ([eb05030](https://github.com/biosimulations/Biosimulations/commit/eb05030d658ff4570eac67c4fb08a91e38d70950))
- **mateiral:** add a material topbar ([4325dcc](https://github.com/biosimulations/Biosimulations/commit/4325dcc4fcb05cfd46a900b2df88bd33ce0ec32f))
- **models:** add query-options model definition ([d628289](https://github.com/biosimulations/Biosimulations/commit/d6282893f47393846af613e1cc4173c0f3b550ef))
- **projects:** add view project ([dab3c59](https://github.com/biosimulations/Biosimulations/commit/dab3c5917c4905c617e1fed78a50854aabc80d37))
- **pwa:** Add pwa capabilities ([27b9050](https://github.com/biosimulations/Biosimulations/commit/27b90508a31e840acc16fbccae809b2b288fde3b))
- **resources:** Resources now have owner embedded ([f03c30b](https://github.com/biosimulations/Biosimulations/commit/f03c30bf7ef7ad5096b27aaf431730c0e402c6ac))
- **serializers:** serializers read files properly ([c0eb5c1](https://github.com/biosimulations/Biosimulations/commit/c0eb5c15dae85aa9e164643e04c957359afca6fa))
- **service:** Breadcrumb service generates breadcrums ([68efc1a](https://github.com/biosimulations/Biosimulations/commit/68efc1a43e0ff52250b16edc73b4351c0b6f647c))
- **services:** add config service ([98ed473](https://github.com/biosimulations/Biosimulations/commit/98ed473bfa04ba40358ef53b376932a660abc286))
- **services:** add file service ([0538686](https://github.com/biosimulations/Biosimulations/commit/053868605725bc0e5fad0944c53502737e46b950))
- **services:** add test organism to metadataservice ([91f2503](https://github.com/biosimulations/Biosimulations/commit/91f25039a9b8d39781bd7b7a2137362de18c7358))
- **services:** Resource services return new model ([0c32ba5](https://github.com/biosimulations/Biosimulations/commit/0c32ba582802824e1dd4e9f64f332e145ead1bee))
- **shared:** add an authentication library ([904cd59](https://github.com/biosimulations/Biosimulations/commit/904cd59b320a33ec8447980e79c12bbcb3121a29))
- **shared:** add fields to remote file ([b85ff9d](https://github.com/biosimulations/Biosimulations/commit/b85ff9d7061a1e62d191b017078e680d5429c7b9))
- **shared:** add more model serializing ([b06129d](https://github.com/biosimulations/Biosimulations/commit/b06129d1b0b7e2b4499c94fb34490b38a15345a3))
- **visualization:** More flexible 2D visualization to better match needs of SED-ML L1V3 and BioModels ([31c96de](https://github.com/biosimulations/Biosimulations/commit/31c96de50f6c1323b574e00f400e5e13ade93060))
- add construction gaurd in prodction mode ([84708a2](https://github.com/biosimulations/Biosimulations/commit/84708a21f656318260e6162518fef84997051916))
- create a debugger component ([fd47c37](https://github.com/biosimulations/Biosimulations/commit/fd47c374f70d4dc28646cec31f995075c778d7a3))
- **shared:** add under constrcution gaurd ([d931f6d](https://github.com/biosimulations/Biosimulations/commit/d931f6de41fc1adfa558355664b651ce4a160c6d))
- **shared:** Remote file has a method to create from File ([f3360cf](https://github.com/biosimulations/Biosimulations/commit/f3360cfb6cfbb733d46501327986e7c9aec565d2))
- **simulations:** Simulations now have embedded models ([4dd512b](https://github.com/biosimulations/Biosimulations/commit/4dd512bc1944a6a2b2e8b00a8d5f48bd61939e9f))
- **users:** Add different snackbars ([056ffe2](https://github.com/biosimulations/Biosimulations/commit/056ffe2e4f26b4f43835a9cb7978e1baa548cf01))
- **visualizations:** add view async capapbility ([d2eb414](https://github.com/biosimulations/Biosimulations/commit/d2eb414a45e9d4eee2caf0edcf87e29c1b066add))

### Reverts

- "Formatted Files. [skip ci]" ([0414aca](https://github.com/biosimulations/Biosimulations/commit/0414aca53ff1c5ee5de9920fdb3e7e5810582a1b))
- Revert "debug redis host" ([7264f87](https://github.com/biosimulations/Biosimulations/commit/7264f87a60c994a655fd80536c33eb234efb5d2b))
- Revert "Feat(Combine-serive): Update API Specification" ([7ada4ce](https://github.com/biosimulations/Biosimulations/commit/7ada4ce33538e4bcd3bae9ee798ea46e417cee99))
- Revert "Updated Ontologies" ([aab6f70](https://github.com/biosimulations/Biosimulations/commit/aab6f708ab3f2d6e39f780bcb4bff1a188376748))
- Revert "Bump @sendgrid/mail from 7.4.1 to 7.4.2 in /biosimulations (#1975)" (#1979) ([a96e0ea](https://github.com/biosimulations/Biosimulations/commit/a96e0ea8f68b97fc2722d464614ee002bcd479c5)), closes [#1975](https://github.com/biosimulations/Biosimulations/issues/1975) [#1979](https://github.com/biosimulations/Biosimulations/issues/1979)
- Revert "Bump @nrwl/cli from 10.4.1 to 11.0.20 in /biosimulations (#1855)" (#1860) ([b7d637e](https://github.com/biosimulations/Biosimulations/commit/b7d637ec1decf056c1642608d821c7b53422b46d)), closes [#1855](https://github.com/biosimulations/Biosimulations/issues/1855) [#1860](https://github.com/biosimulations/Biosimulations/issues/1860)
- Revert "merging enumerations of simulation status; aligning names of properties 'resultsSize' and 'resultSize'" ([08e96f2](https://github.com/biosimulations/Biosimulations/commit/08e96f2fe8c7e88c7fbb2a6653416182fc438700))
- Revert "Revert "styling nested lists and lists after paragraphs"" ([16a2d5b](https://github.com/biosimulations/Biosimulations/commit/16a2d5ba5e39a79c64bca5af77c1382c76dfb2c5))
- Revert "Revert "adding management for app-specific configuration (e.g,, appName, logo, etc."" ([1c6190a](https://github.com/biosimulations/Biosimulations/commit/1c6190ab50c95359233e04befc2c778fa1dbd5c1))
- Revert "Revert "adding documentation"" ([f71cfcd](https://github.com/biosimulations/Biosimulations/commit/f71cfcddcf15006f3f84660517c7fdf52ada0dbe))
- Revert "Revert "editing help"" ([4aea612](https://github.com/biosimulations/Biosimulations/commit/4aea6128b0f4fbb2d881efbcc5eb2416527f05c9))
- Revert "Revert "adding documentation of supported SED-ML features"" ([32beb34](https://github.com/biosimulations/Biosimulations/commit/32beb34d5a823e3975d6115d8f225044fe6297a5))
- Revert "changing completed to updated" ([d9fe18e](https://github.com/biosimulations/Biosimulations/commit/d9fe18edbbf4a1d0307caecdeb7a6927b0338abf))
- Revert "update to angular 10" ([e1b9fcf](https://github.com/biosimulations/Biosimulations/commit/e1b9fcfae21f5bd8b19e4cee23a8df2b82c1df02))
- Revert "Bump husky from 4.0.1 to 4.0.5 in /CRBM-Viz (#290)" ([4dd2cc3](https://github.com/biosimulations/Biosimulations/commit/4dd2cc3eee2fcf33bb24ca8ef999162823e875a7)), closes [#290](https://github.com/biosimulations/Biosimulations/issues/290)

- feat (dispatch-api): remove download endpoint ([785ad27](https://github.com/biosimulations/Biosimulations/commit/785ad27f1477ac6122c2a735cb46201928a0f754))
- feat (dispatch) : Change datamodel of returned results ([ba42dcc](https://github.com/biosimulations/Biosimulations/commit/ba42dcc8b74068a280c6dc2f7d915c8c27a55f45))

### BREAKING CHANGES

- The /download endpoint has been removed. Should be replaced by /results/download
- The results are now returned as an array of objects (AOS) rather than an object of arrays (SOA)
# Changelog

## [8.8.0](https://github.com/biosimulations/biosimulations/compare/v8.7.1...v8.8.0) (2022-01-08)


### Features

* added workflow to delete temporary COMBINE archives ([0d722e9](https://github.com/biosimulations/biosimulations/commit/0d722e939b19726c0f91c310d8a0a6a31217aed8))
* **api,combine-api,dispatch,platform:** added support for references for projects ([a544969](https://github.com/biosimulations/biosimulations/commit/a5449691805801726b643285a018dbff77e06a00))

## [8.7.1](https://github.com/biosimulations/biosimulations/compare/v8.7.0...v8.7.1) (2022-01-06)


### Bug Fixes

* add script ignore false to sharp install ([9b90d9b](https://github.com/biosimulations/biosimulations/commit/9b90d9bfeefe1e8b1044e9e2685e5de15bb3b2e4))
* **dispatch-service:** added dependency for sharp to Dockerfile ([ba344c7](https://github.com/biosimulations/biosimulations/commit/ba344c7649c97b842b9794e3f58bd8af3cbe2712))
* **dispatch,platform:** fixed name and URL for log format in files tab ([f13647f](https://github.com/biosimulations/biosimulations/commit/f13647fcaac05b5abdfb5ad4a754c3f4c18586ac))

## [8.7.0](https://github.com/biosimulations/biosimulations/compare/v8.6.0...v8.7.0) (2022-01-06)


### Bug Fixes

* **dispatch-service:** separated input and output files for simulation runs ([dc98a46](https://github.com/biosimulations/biosimulations/commit/dc98a4604add2541e1214afa4701eb829c634a3f))
* **dispatch,platform,ui:** corrected layout of metadata columns ([8044bc6](https://github.com/biosimulations/biosimulations/commit/8044bc62d06057696b1a12d22d61a9287d1dbd7f))
* **dispatch,platform:** added modeling methods to metadata ([6da720e](https://github.com/biosimulations/biosimulations/commit/6da720ef4504515f9d9eafdb11e8a76d56a448cc))
* **dispatch:** add rel noopenor for external links ([9efa059](https://github.com/biosimulations/biosimulations/commit/9efa059a50435583859f17ff05476694a84a5e59))
* **ui:** fixed setting of open control panel in table controls ([ca529a3](https://github.com/biosimulations/biosimulations/commit/ca529a3be58ca72ab4791732de69ff0efa8383fa))


### Features

* **combine-api:** added utility methods for reading S3 files ([4e3df4f](https://github.com/biosimulations/biosimulations/commit/4e3df4fee4247f1d4ae371f634cb6af8084b93ee))
* **dispatch-service:** added support for SLURM constraints ([df2acba](https://github.com/biosimulations/biosimulations/commit/df2acba2839c4def3139dc0a889af4242889fbde))
* **dispatch,platform,ui:** added markdown rendering for project descriptions ([818d267](https://github.com/biosimulations/biosimulations/commit/818d2673f20c1a87f235618cae242c129b49630f))
* **dispatch,platform,ui:** interleaved metadata about files into files tab ([9807406](https://github.com/biosimulations/biosimulations/commit/9807406cbcfc1c78b793d2068a98b4693725741e))
* **dispatch:** made it easier to get errors with simulation projects ([4c7635f](https://github.com/biosimulations/biosimulations/commit/4c7635fb9997d8cb33f4c2051bb6ffeb9958c42e))
* **ontology:** added additional formats used by Physiome ([3e951d0](https://github.com/biosimulations/biosimulations/commit/3e951d02cd6633aca35d444ce60dd076dbcb1dc8))


### Performance Improvements

* **api:** added caching for getting properties of ontology terms ([c273e88](https://github.com/biosimulations/biosimulations/commit/c273e88cb31efeaf7872aa48ad1aa38a05524b80))
* **dispatch,platform:** reduced thumbnail image sizes ([049af36](https://github.com/biosimulations/biosimulations/commit/049af36affaf05b56a35834c0c342803d74ba46c))

## [8.6.0](https://github.com/biosimulations/biosimulations/compare/v8.5.6...v8.6.0) (2021-12-31)


### Features

* **ontology:** added formats used by Physiome model repository ([b407218](https://github.com/biosimulations/biosimulations/commit/b40721893c8592e1418f21fe78462d1978270056))

## [8.5.6](https://github.com/biosimulations/biosimulations/compare/v8.5.5...v8.5.6) (2021-12-27)


### Bug Fixes

* update package lock version ([b91c669](https://github.com/biosimulations/biosimulations/commit/b91c669745abee0599cfd002fd18f557acd8e155))

## [8.5.5](https://github.com/biosimulations/biosimulations/compare/v8.5.4...v8.5.5) (2021-12-23)


### Bug Fixes

* **api:** tried to correct logging for errors in uploading archives ([b492f02](https://github.com/biosimulations/biosimulations/commit/b492f02c04da3bc6bf1bdd7ec022ce85187ddfa1))

## [8.5.4](https://github.com/biosimulations/biosimulations/compare/v8.5.3...v8.5.4) (2021-12-23)


### Bug Fixes

* improved error logging, increasing AWS timeout ([cccf802](https://github.com/biosimulations/biosimulations/commit/cccf802c76f154ed09dc53a7b1a41df286d935a8))

## [8.5.3](https://github.com/biosimulations/biosimulations/compare/v8.5.2...v8.5.3) (2021-12-23)


### Bug Fixes

* **dispatch-service:** fixed retrying of files and specs to avoid conflicts on incomplete posts ([620f0eb](https://github.com/biosimulations/biosimulations/commit/620f0eb59a1f347ea9bbbaaf0b3d0276654d8611))

## [8.5.2](https://github.com/biosimulations/biosimulations/compare/v8.5.1...v8.5.2) (2021-12-23)


### Bug Fixes

* **api,dispatch-service:** fixed logging for complete processor failures ([f2119d4](https://github.com/biosimulations/biosimulations/commit/f2119d46741eb85f89a91d84eeee93903900acc8))

## [8.5.1](https://github.com/biosimulations/biosimulations/compare/v8.5.0...v8.5.1) (2021-12-23)


### Bug Fixes

* **dispatch-service:** retry project publication after run completion ([e686664](https://github.com/biosimulations/biosimulations/commit/e68666451fc61bd27dd25e3775a4477fca5385e1))

## [8.5.0](https://github.com/biosimulations/biosimulations/compare/v8.4.1...v8.5.0) (2021-12-22)


### Bug Fixes

* **simulators-api:** fixed updating of updated timestamp; addresses [#3878](https://github.com/biosimulations/biosimulations/issues/3878) ([13ee425](https://github.com/biosimulations/biosimulations/commit/13ee42542af3f15c13954e7873313581cd76647c))
* **ui,dispatch:** fixed unselecting files; closes [#3875](https://github.com/biosimulations/biosimulations/issues/3875) ([2f84e18](https://github.com/biosimulations/biosimulations/commit/2f84e187ee1ea2c865925d486b5235682814e0da))


### Features

* **ui:** added instructions to refresh on failures ([0181c3b](https://github.com/biosimulations/biosimulations/commit/0181c3be83cfe5602d258ca8276599348c123f9e))

## [8.4.1](https://github.com/biosimulations/biosimulations/compare/v8.4.0...v8.4.1) (2021-12-22)


### Bug Fixes

* **combine-api:** fixed reading shared configuration ([5a83e6e](https://github.com/biosimulations/biosimulations/commit/5a83e6e0f4802eed023781dac8881f4c8e251a2e))

## [8.4.0](https://github.com/biosimulations/biosimulations/compare/v8.3.0...v8.4.0) (2021-12-22)


### Bug Fixes

* added project id, owner to CompleteJob for failures ([0a8d834](https://github.com/biosimulations/biosimulations/commit/0a8d834e100fc396da08403e5cef6e8c703c133c))
* **api,dispatch:** fixed data model for simulation results ([7cd2b5e](https://github.com/biosimulations/biosimulations/commit/7cd2b5e2a3d37e7544c16b383ee6594b65afdf69))
* **api,dispatch:** fixed data type for simulation results ([610cbc9](https://github.com/biosimulations/biosimulations/commit/610cbc97530b8b9965a98a356beb23e27fbf3cdc))
* **api:** changed file URL validation to allow un-encoded URLs ([f78bc43](https://github.com/biosimulations/biosimulations/commit/f78bc43fd89bb990841786eacb5118b5645628c4))
* fixed external simulators API endpoint ([25b3a3c](https://github.com/biosimulations/biosimulations/commit/25b3a3c0aabbf14949cf5386438b4bee87c95a03))
* fixed spinner for loading table data ([e2e1314](https://github.com/biosimulations/biosimulations/commit/e2e13149988f7cd74f95849f9c473c914d44f5ee))


### Features

* added checks that S3 files were deleted ([390300c](https://github.com/biosimulations/biosimulations/commit/390300c6039cede0176f7a5824b01687dbae5a94))
* **api:** added cache for project summaries ([6fc5bb7](https://github.com/biosimulations/biosimulations/commit/6fc5bb705a3cadb3ee2379145b3a362f8abac354))
* **api:** added checks that S3 files were deleted ([5d3ccb9](https://github.com/biosimulations/biosimulations/commit/5d3ccb91f964b11e9ffca1db1cb5f51d2c1d4389))
* **api:** added handled for NaN and Inf from HSDS ([7802310](https://github.com/biosimulations/biosimulations/commit/78023109c8a0b3c1c6d4e80ec65aec680b8a6172))
* **combine-api:** relaxed required metadata for simulation projects ([2be28f1](https://github.com/biosimulations/biosimulations/commit/2be28f1a896256e29ea7f9387e66f6283538e094))
* directed href targets ([d609f80](https://github.com/biosimulations/biosimulations/commit/d609f80c71549213eaac6b9b1fd8a43b7d2c038c))

## [8.3.0](https://github.com/biosimulations/biosimulations/compare/v8.2.1...v8.3.0) (2021-12-20)


### Bug Fixes

* **dispatch-service:** fix dispatch service post limit ([95653f8](https://github.com/biosimulations/biosimulations/commit/95653f88710f37babc19ccfa85b7c7caabc9fec3)), closes [#3828](https://github.com/biosimulations/biosimulations/issues/3828)
* **dispatch:** correct security issue with untrusted html input ([8b464b6](https://github.com/biosimulations/biosimulations/commit/8b464b64ecc1992031e49b76a4765ce68a70b7f4))
* **ui:** fixed spinner exit for table component ([f26802c](https://github.com/biosimulations/biosimulations/commit/f26802cf284df664509f2052a5f29a4443ee9a36))


### Features

* **api:** added project summary caching at creation and updating ([20d99f2](https://github.com/biosimulations/biosimulations/commit/20d99f2bcccba153a4aae10b9872c212935563f6))
* **dispatch:** added example simulation runs for Brian 2 ([710df68](https://github.com/biosimulations/biosimulations/commit/710df689dc79f3ac837ddda2034398371eb5f084))

## [8.2.1](https://github.com/biosimulations/biosimulations/compare/v8.2.0...v8.2.1) (2021-12-16)


### Bug Fixes

* add gtag snippet to dispatch and simulators ([f1b6332](https://github.com/biosimulations/biosimulations/commit/f1b633298b064a3a9d6ce8bdc404fd815ead5a5c))
* **config:** add default server limit to config ([496b430](https://github.com/biosimulations/biosimulations/commit/496b430efd96e7d2b13b102ea8f7ef9d25b8e35a)), closes [#3828](https://github.com/biosimulations/biosimulations/issues/3828)

## [8.2.0](https://github.com/biosimulations/biosimulations/compare/v8.1.0...v8.2.0) (2021-12-16)


### Bug Fixes

* fixed log validation ([86fc30d](https://github.com/biosimulations/biosimulations/commit/86fc30d388cf1ee170952be7efabdbd5bc5faca7))


### Features

* add angular analytics package ([3363f7b](https://github.com/biosimulations/biosimulations/commit/3363f7bcae10fcc07383e36c617bb960b054f380))
* added implementation of analytics and user consent ([2d87bb1](https://github.com/biosimulations/biosimulations/commit/2d87bb16abe21c9af02f402e3cfdb9265b3605e6))
* **dispatch,platform,simulators:** add cookie consent and privacy settings to frontend apps ([e84cdea](https://github.com/biosimulations/biosimulations/commit/e84cdeaa8a230a068afbb490dc33a796a441cc59))

## [8.1.0](https://github.com/biosimulations/biosimulations/compare/v8.0.0...v8.1.0) (2021-12-16)


### Bug Fixes

* fixed display of files in subdirectories ([ee62bbe](https://github.com/biosimulations/biosimulations/commit/ee62bbebf4cdeafcc3ec24f20a16bbc1178a0a82))
* fixed file size extraction for empty files ([b3dce39](https://github.com/biosimulations/biosimulations/commit/b3dce39fc83d9d846d1e1a87e5435e9cdc5078b2))
* fixed simulators view endpoint method ([e40f4cd](https://github.com/biosimulations/biosimulations/commit/e40f4cd29497f096d7341fd3c7b9c589b27f0093))
* set minimum time step to 1 ([e92ffa5](https://github.com/biosimulations/biosimulations/commit/e92ffa5de8a163b7e844a8f47cdd55e02b2b3155))


### Features

* **combine-api:** added error messages for invalid S3 bucket configuration ([626a2f2](https://github.com/biosimulations/biosimulations/commit/626a2f2a6deadaaa71517a09c19451f975316ad0))

## [8.0.0](https://github.com/biosimulations/biosimulations/compare/v7.0.0...v8.0.0) (2021-12-15)


### Bug Fixes

* **config:** change name of env variable to avoid clash ([caf69f4](https://github.com/biosimulations/biosimulations/commit/caf69f4ef11439ad6b6269ee41e904283219225c))
* **config:** correct endpoints for s3 contents path. Add some testing ([73f6034](https://github.com/biosimulations/biosimulations/commit/73f6034561fd88690a89f45d76dde9ed681b55ee)), closes [#3755](https://github.com/biosimulations/biosimulations/issues/3755)
* **config:** fix endpoints for ontology url ([f7ba9c5](https://github.com/biosimulations/biosimulations/commit/f7ba9c5f55048092280a0eeef0740bccc84d4c4e)), closes [#3771](https://github.com/biosimulations/biosimulations/issues/3771)
* **config:** fix prod file url ([5211afd](https://github.com/biosimulations/biosimulations/commit/5211afd21258dedd1483b7509270accfdc6f8dc8))
* **config:** fix storage health endpoint ([14d303c](https://github.com/biosimulations/biosimulations/commit/14d303c8c304468238b446acb1cbb7a4023e29b3))
* **dispatch-service:** corrected storage endpoint in sbatch sevice to external ([8eb34d4](https://github.com/biosimulations/biosimulations/commit/8eb34d4091e010b8c7b4ada5d29ceb1bf5c5a4d5))
* fixed broken links in documentation ([fe96ded](https://github.com/biosimulations/biosimulations/commit/fe96dedddd943974a5768a0f4439095f6bea8958))
* fixed broken links in documentation ([8de8a6f](https://github.com/biosimulations/biosimulations/commit/8de8a6f40ae046a95d5557d7d453220037eb9e27))
* removed invalid and unecessary workflow_call secret inputs ([b925dfe](https://github.com/biosimulations/biosimulations/commit/b925dfed75d24546ac5a01c6ac43fbf378b594b7))
* **simulators-api:** fix permissions for deletion of all simulation runs ([24aef29](https://github.com/biosimulations/biosimulations/commit/24aef2972670e4046707fb98072b0834f8435472)), closes [#3767](https://github.com/biosimulations/biosimulations/issues/3767)


### Code Refactoring

* **api:** simplify management of files ([e16c35f](https://github.com/biosimulations/biosimulations/commit/e16c35f374f5cbf0acda5f26288a1ed5c1ce04ec))


### Features

* **api,dispatch-service,dispatch:** expanded to full SED-ML data model ([2550def](https://github.com/biosimulations/biosimulations/commit/2550defc9918ae44dd2d5df53b6e2035ebcb7a00))
* **api:** improved error logging ([4d75193](https://github.com/biosimulations/biosimulations/commit/4d75193c68578f4dd48740795dfe9a630109366e))
* **auth:** add scope for deleting all simulation runs ([77a859c](https://github.com/biosimulations/biosimulations/commit/77a859caec20639adf77e6b0ac6d5ad00ecbb1ca))
* **combine-api:** expanded to full Python SED-ML data model ([912f8d5](https://github.com/biosimulations/biosimulations/commit/912f8d57bf205c76f7a12b4147dbf63e7bed0e89))
* **combine-api:** updated dependencies ([a7bf6e9](https://github.com/biosimulations/biosimulations/commit/a7bf6e9e4060658f62837e74172169abeb217f80))
* **config:** load endpoints dynamically if not in browser ([90ddc18](https://github.com/biosimulations/biosimulations/commit/90ddc1895fb71d04663c86965bf8fa3a6d1ec785)), closes [#3585](https://github.com/biosimulations/biosimulations/issues/3585)
* **dispatch-service:** added env variables to enable simulators to get number of CPUs ([fbdf8fa](https://github.com/biosimulations/biosimulations/commit/fbdf8fa354a6b06339d9831385c649a0e4d02e24))
* **dispatch,dispatch-service:** added support for multidimensional plots ([6681d25](https://github.com/biosimulations/biosimulations/commit/6681d25aed91b2bc51f5442ab929390765c9ed63))
* **dispatch,platform,api:** added filter and display of SED-ML file provenance ([d2a6c29](https://github.com/biosimulations/biosimulations/commit/d2a6c29a77ea1fc10289e576f10a6ef3bbb467b9))
* **dispatch:** extend vega export to multidimensional data ([4474f0c](https://github.com/biosimulations/biosimulations/commit/4474f0cf3cc7711b5e7b5dd9062e8d8e1c2585d9))
* expanded COMBINE archive creation to all types of model changes ([f65ffc5](https://github.com/biosimulations/biosimulations/commit/f65ffc5e6dde5a1edbbebbee83a20a73b1305479))


### Performance Improvements

* **api:** extract the combine archive directly on s3 via streams ([b5a0f08](https://github.com/biosimulations/biosimulations/commit/b5a0f0869610dc73e5119ecf44054e35f8354e15)), closes [#3094](https://github.com/biosimulations/biosimulations/issues/3094)


### Reverts

* "chore(deps): update dependency typescript to v4.5.4" ([f40932e](https://github.com/biosimulations/biosimulations/commit/f40932e6a41f56d94343ccad0af54d7368012572))
* "chore(deps): update typescript-eslint monorepo to v5.7.0" ([dee846f](https://github.com/biosimulations/biosimulations/commit/dee846f58f5989ef2f40733d97cbbff2b5f0722c))
* "chore(deps): update typescript-eslint monorepo to v5.7.0" ([16bcd80](https://github.com/biosimulations/biosimulations/commit/16bcd80c16b51ca2b3471b62f3fe10e0e2b448c7))


### BREAKING CHANGES

* **api:** the download project endpoint will fail for all previously submitted simulation
runs. Runs submitted prior to this change will not be retrievable by the api or applications

## [7.0.0](https://github.com/biosimulations/biosimulations/compare/v6.1.0...v7.0.0) (2021-12-02)


### Bug Fixes

* **api,dispatch-service,hsds:** fixed retrying in APIs ([b3109a1](https://github.com/biosimulations/biosimulations/commit/b3109a1fa5d37d0a65134adde6d9aaef90b19ac0))
* **config:** reverted changes to localhost ([04f2c42](https://github.com/biosimulations/biosimulations/commit/04f2c42263b8b98facda48a712bfb9f7e2d2e50d))
* corrected BioSimulators auth audience ([5c8b18d](https://github.com/biosimulations/biosimulations/commit/5c8b18d0d973bb39c8ba34375a9ba33219a76b44))
* corrected BioSimulators auth audience ([fce4347](https://github.com/biosimulations/biosimulations/commit/fce434727d36d28f048ea3b2d852ab56574a6768))
* corrected filtering for numerical columns ([d4bf35f](https://github.com/biosimulations/biosimulations/commit/d4bf35f09605fe5210c5b3ae953e0f6156624bec))
* **dispatch:** fix displaying of data in cases where simulation fails and metadata is not present ([da18704](https://github.com/biosimulations/biosimulations/commit/da18704a8b8fef837b9a835d12e95c6f6644d237)), closes [#3705](https://github.com/biosimulations/biosimulations/issues/3705)
* **dispatch:** fix file uploading ([30f9cb9](https://github.com/biosimulations/biosimulations/commit/30f9cb93e6ed9c52b08a62655de5693bd9552325)), closes [#3719](https://github.com/biosimulations/biosimulations/issues/3719)
* **dispatch:** fix uploading of files ([f597e90](https://github.com/biosimulations/biosimulations/commit/f597e9050b7385b8e17b7179476512d15cb3d723)), closes [#3719](https://github.com/biosimulations/biosimulations/issues/3719)
* fixed validation by upgrading from broken version of class-transformer ([ca20307](https://github.com/biosimulations/biosimulations/commit/ca203077941683837fb5426d1a7eb619d49d6b4c))
* **platform,ui:** fixed project browse for mobile ([fcfb552](https://github.com/biosimulations/biosimulations/commit/fcfb552a52b641d9f616f1e0785929d1b4bf413e))
* **platform:** fixed position of seach/filter button ([133f08c](https://github.com/biosimulations/biosimulations/commit/133f08c99c691b4cb62ba5c548abb2d7a41a7013))
* **ui:** fixed autocomplete filter ([cb26564](https://github.com/biosimulations/biosimulations/commit/cb26564c503e9d10f4423c21e27993e7c77e1802))


### Features

* added ability to filter and search projects ([21c14b4](https://github.com/biosimulations/biosimulations/commit/21c14b4fec4c1a0e56773b16e7b9e9e46c306636))
* added new component to enable components for routes to push buttons into the breadcrumbs area ([7918db1](https://github.com/biosimulations/biosimulations/commit/7918db1fccba56f3b9127bccdf042ddb7e7d83c4))
* **api,dispatch,dispatch-service:** expanded support for failed simulation runs ([9e7e71c](https://github.com/biosimulations/biosimulations/commit/9e7e71c80a0a81de3195648304456dfcf293b00c))
* **api,platform:** added owners, organizations to project view with hyperlinks ([63d9457](https://github.com/biosimulations/biosimulations/commit/63d9457a5244212523b727c3f7a01117b7253711))
* **api,platform:** started to display ownership of projects ([6f8378a](https://github.com/biosimulations/biosimulations/commit/6f8378a14d2f3d3291b36d11ada2581d1e6bd2be))
* **api:** add check for data service to status check ([d8fbbc5](https://github.com/biosimulations/biosimulations/commit/d8fbbc5308915d23e23dcd9b5962c869236f8988)), closes [#3649](https://github.com/biosimulations/biosimulations/issues/3649)
* **api:** added new scope for externally validating simulation runs ([e3bd698](https://github.com/biosimulations/biosimulations/commit/e3bd698c989c48cda49e242d2677859d928b9324))
* **api:** added URLs to accounts, organizations ([c6926b9](https://github.com/biosimulations/biosimulations/commit/c6926b95abbd8873bc4b56fbe77038b2cc2a76be))
* **api:** began to limit publication requests to model repositories ([e830002](https://github.com/biosimulations/biosimulations/commit/e8300028c8f8e22a2ef59938b7c83880167df5a7))
* **api:** working on restricting requests for publication with simulation run requests ([af7b6a0](https://github.com/biosimulations/biosimulations/commit/af7b6a0d46cedff4741d07502ccf8ef07fbbaa89))
* **dispatch:** clarified units of columns of simulation runs table ([3ae498e](https://github.com/biosimulations/biosimulations/commit/3ae498ecc1902c6fd688889427b0fdf28ad25b2e))
* improved table searching for data with accents ([d33617b](https://github.com/biosimulations/biosimulations/commit/d33617b73d1968f74cacaeae0bdd92e2a3aa3339))
* **platform:** added filter for publication status ([509eaa4](https://github.com/biosimulations/biosimulations/commit/509eaa4a431eb481cde071edc666ee7bdbe3961f))
* **platform:** scroll to top on opening projects search/filter ([26de2b6](https://github.com/biosimulations/biosimulations/commit/26de2b6f1470cc9fbb661447ac0542c2c7ea3483))
* **platform:** started to add filtering and searching for projects ([2505d7f](https://github.com/biosimulations/biosimulations/commit/2505d7f2bdca62c3b24d4b7252cb2513a9b5aa1b))
* **simulators:** expanded simulators filters ([1a76b7e](https://github.com/biosimulations/biosimulations/commit/1a76b7e5934c30793673ac8b2250c5b370ba2fa9))
* **ui,platform:** added autocomplete filter for attributes with many values ([049869f](https://github.com/biosimulations/biosimulations/commit/049869f919dcf9b714fc2d27389d1b1b6ce21971))
* **ui:** add custom caruousel component ([9a61d4f](https://github.com/biosimulations/biosimulations/commit/9a61d4f45592aa0b4c6d7205171a7b3a87b0330f))
* **ui:** replace npn-slider with custom component ([1cc79d5](https://github.com/biosimulations/biosimulations/commit/1cc79d54bc20114247ee6af36e2d3ade44c2f64d)), closes [#3706](https://github.com/biosimulations/biosimulations/issues/3706)


### Performance Improvements

* **api,dispatch,dispatch-service,simulators-api:** removed unnecessary return of new resources ([5fce07f](https://github.com/biosimulations/biosimulations/commit/5fce07f6b52385ab5d04c276b635f78aa06e2eea))
* **api:** don't return the log after creating ([5e37c6e](https://github.com/biosimulations/biosimulations/commit/5e37c6eca3aea9420467884f75f05558b6d30c0b)), closes [#3609](https://github.com/biosimulations/biosimulations/issues/3609)


### Reverts

* Revert "refactor(auth): removed auth/open endpoint" ([f983628](https://github.com/biosimulations/biosimulations/commit/f983628bd2b048b2a9e2b3875c4508f7a7b62746))
* "refactor: cleaned up building Angular apps" ([d21b2ed](https://github.com/biosimulations/biosimulations/commit/d21b2edc4530ff644d92ac6f71385dc83052fa3d))
* revert "refactor: organized endpoints configuration" ([a5e93c3](https://github.com/biosimulations/biosimulations/commit/a5e93c38a268995474525dc71e1fb72ef8fdf968)), closes [#3625](https://github.com/biosimulations/biosimulations/issues/3625)
* revert change to build front end apps ([5870fbf](https://github.com/biosimulations/biosimulations/commit/5870fbf953b95bd55d6def44aac7a5628c0a7265))


### BREAKING CHANGES

* **api:** The logs post endpoint no longer returns the log that was created. For that, use a
GET request after posting the log.

## [6.1.0](https://github.com/biosimulations/biosimulations/compare/v6.0.2...v6.1.0) (2021-11-14)


### Bug Fixes

* **auth:** fix import ([a7ba6b8](https://github.com/biosimulations/biosimulations/commit/a7ba6b87404b621d22c90c9fb37164b836d39f85))
* **auth:** handle case of no custom permissions ([1c3d760](https://github.com/biosimulations/biosimulations/commit/1c3d760f42ea23bdb3c26a42dce81b09189f7f92))
* corrected capitalization of BioSimulations ([3d981ee](https://github.com/biosimulations/biosimulations/commit/3d981ee737b852488dfd7ff8aba1317c11b7f236))
* debugged testing COMBINE API ([35672ce](https://github.com/biosimulations/biosimulations/commit/35672ceb0ce44474f8644c1bef71584208d778c7))
* debugged testing COMBINE API ([d18cb86](https://github.com/biosimulations/biosimulations/commit/d18cb864f29500229fa24666d9525bc228324476))
* debugged testing COMBINE API ([2a1e6b3](https://github.com/biosimulations/biosimulations/commit/2a1e6b3fc24fcdba88f7cc839d3655d6b4be80ce))
* debugged testing COMBINE API ([c7e1cb9](https://github.com/biosimulations/biosimulations/commit/c7e1cb9efc3daacaa6eddc1fd48806d25ec9927b))
* **dispatch:** corrected run URLs in check simulation run tool ([a1894fa](https://github.com/biosimulations/biosimulations/commit/a1894fae9b406a62fb4d10a30d01c9988d76a95a))
* **dispatch:** fixed dispatch simulation run view; closes [#3088](https://github.com/biosimulations/biosimulations/issues/3088) ([26a8d5a](https://github.com/biosimulations/biosimulations/commit/26a8d5a49b577c89eb04cd4b58314122f2d844c1))
* **dispatch:** fixed highlight.js import for log formatting ([0320cc2](https://github.com/biosimulations/biosimulations/commit/0320cc201838549a4d786c117571db132ddad600))
* fixed links, warnings ([ecb68fe](https://github.com/biosimulations/biosimulations/commit/ecb68febccdb4f17308aa62f0172adfcc7c68554))
* fixed python code highlighting ([ef88e12](https://github.com/biosimulations/biosimulations/commit/ef88e1248bbf2f9c8493c6f9a94a7c359d4497d0))
* fixed typos, added spelling exceptions ([7c77dc2](https://github.com/biosimulations/biosimulations/commit/7c77dc2db2b38563fea78118b1b65adca0510198))
* removed example with COMBINE archive that intentionally fails ([04b4c6e](https://github.com/biosimulations/biosimulations/commit/04b4c6ed1923ad2fa8ec1ae1f963e638cbe95807))


### Features

* **api,dispatch-service:** improved error messages and retrying ([a6b2693](https://github.com/biosimulations/biosimulations/commit/a6b26935a12e66a5bc3ab0b00713d8613c7a3f5d))
* **api:** improved reporting of errors with inconsistent data ([1c6c223](https://github.com/biosimulations/biosimulations/commit/1c6c2233c353f534efa0503f63b314adb9b97738))
* **dispatch-service:** add logging to processing posts to api ([1ca1900](https://github.com/biosimulations/biosimulations/commit/1ca19009cd33f6fbbb79ba9fb15779e09e20c1c8))
* **dispatch-service:** add retries for posting processing results ([fe9cddc](https://github.com/biosimulations/biosimulations/commit/fe9cddce4c846f15fdfd002bb8465e0e48bb3bd8)), closes [#3531](https://github.com/biosimulations/biosimulations/issues/3531)
* improved docs ([ab6722b](https://github.com/biosimulations/biosimulations/commit/ab6722b558a0b4e0d5aca265c84a1d0afd9f2558))
* improved docs ([209a421](https://github.com/biosimulations/biosimulations/commit/209a421c5c177600dd3e88c9469020db3c2fa51c))

## [6.0.2](https://github.com/biosimulations/biosimulations/compare/v6.0.1...v6.0.2) (2021-11-10)


### Bug Fixes

* **hsds:** update the hsds client ([419bfe9](https://github.com/biosimulations/biosimulations/commit/419bfe9a7ef4b81eea9400e2a4aa7587d734b4bf)), closes [#3317](https://github.com/biosimulations/biosimulations/issues/3317)

## [6.0.1](https://github.com/biosimulations/biosimulations/compare/v6.0.0...v6.0.1) (2021-11-10)


### Bug Fixes

* **api:** fixed IsImageDigest validator for non-strings ([53501b1](https://github.com/biosimulations/biosimulations/commit/53501b1eb3f323658aa23ae9007d625f212f6c6f))

## [6.0.0](https://github.com/biosimulations/biosimulations/compare/v5.9.0...v6.0.0) (2021-11-04)


### Bug Fixes

* **api:** build fix for new axios types ([e7ea984](https://github.com/biosimulations/biosimulations/commit/e7ea9849f72eefb1a051c316601e9463c3f74b1b))
* **dispatch-service:** add a temporary check for mistructured logs ([a67bfa1](https://github.com/biosimulations/biosimulations/commit/a67bfa1347b827cc2162080b0167412e33afa12f)), closes [#3482](https://github.com/biosimulations/biosimulations/issues/3482) [#3482](https://github.com/biosimulations/biosimulations/issues/3482)
* **dispatch-service:** remove ssl skip when downloading archive ([b4bd2c0](https://github.com/biosimulations/biosimulations/commit/b4bd2c0553bfa253002efb4063f3be30b321c15c)), closes [#3092](https://github.com/biosimulations/biosimulations/issues/3092)
* **dispatch,api,dispatch-service:** fixed data model for exceptions in simulation run logs ([191f1d3](https://github.com/biosimulations/biosimulations/commit/191f1d36ffc8824eb900b4571b937f0c9191b258))
* fix imports ([a076005](https://github.com/biosimulations/biosimulations/commit/a076005a50f89fa1f5cc76d66ededc44d72633ae))
* **simulators:** fixed text overflow of simulator test results ([7a8b11a](https://github.com/biosimulations/biosimulations/commit/7a8b11a85fd8a5975e39763cfb9247801ac8f5d6))


### Code Refactoring

* **api:** cleaned up simulation run files ([0213fb5](https://github.com/biosimulations/biosimulations/commit/0213fb508b66665b98c64c713ab146505695e2b9))


### Features

* added endpoints for getting summaries of projects ([2e9a54e](https://github.com/biosimulations/biosimulations/commit/2e9a54ebba20f9883bbf2c3dd6297db906da1285))
* **api:** added database model and validation for logs ([6a3a344](https://github.com/biosimulations/biosimulations/commit/6a3a344aa14fbd4952758f676160a9ba472b071e))
* **api:** added endpoints for getting individual SED elements; closes [#3439](https://github.com/biosimulations/biosimulations/issues/3439) ([8aa64fc](https://github.com/biosimulations/biosimulations/commit/8aa64fc205494d27c4f444c9d669031dd3ab49cc))
* **combine-api:** add dynamic module for combine api-client ([e7c2448](https://github.com/biosimulations/biosimulations/commit/e7c24486db6f7bf5b4cbc2b9e9c2c6a67a34bc1b)), closes [#3180](https://github.com/biosimulations/biosimulations/issues/3180)
* **datamodel:** add validation for image digests ([6106e54](https://github.com/biosimulations/biosimulations/commit/6106e54f5d54792b24fc5994a7ac1a42bed0c790))
* **dispatch-service:** enhanced tracking of processing results ([b4f01e3](https://github.com/biosimulations/biosimulations/commit/b4f01e3f4428b64ddb8aa3671fe738ff024515cc))
* **dispatch:** updated publication form, finished switching to Endpoints ([a5bef72](https://github.com/biosimulations/biosimulations/commit/a5bef7248174dec74154087b92bcc6005fa87726))
* **exceptions:** improve error handling ([4a3e8c7](https://github.com/biosimulations/biosimulations/commit/4a3e8c78dcb35756172da30ef803d2954f264bc6))
* **hsds:** handle transient hsds query failures ([f4a19f5](https://github.com/biosimulations/biosimulations/commit/f4a19f531c80ec50d324dec009132450dd74edb2)), closes [#3413](https://github.com/biosimulations/biosimulations/issues/3413)
* **ontology:** added parent/child relationships to ontology terms ([4107f1d](https://github.com/biosimulations/biosimulations/commit/4107f1d3d1f17a7509262caae998fd0222ebf0c3))
* **simulators-api:** add validation for api models ([ce2c5bb](https://github.com/biosimulations/biosimulations/commit/ce2c5bb0616600289e94dc16a7dcc17cfcb27dc4))
* **simulators,dispatch,platform:** added status bar to bottom of apps; closes [#3210](https://github.com/biosimulations/biosimulations/issues/3210) ([3630c23](https://github.com/biosimulations/biosimulations/commit/3630c23ef11325b4bb47012e4ca58ec7fefb6b7c))


### Reverts

* **dispatch-service:** revert using new config service to provide basepath ([e479d2e](https://github.com/biosimulations/biosimulations/commit/e479d2e531383e55ed64823b98b41d7be9130f7f))


### BREAKING CHANGES

* **api:** moves simulation run file information from 'Simulation Files' collection

## [5.9.0](https://github.com/biosimulations/biosimulations/compare/v5.8.0...v5.9.0) (2021-10-24)


### Bug Fixes

* **deps:** update to nx v13 ([48dbf7c](https://github.com/biosimulations/biosimulations/commit/48dbf7cfe4002aed9fcc06237f9bc995539573c2))
* **exceptions:** handle cases of payload too large errors ([522cec8](https://github.com/biosimulations/biosimulations/commit/522cec85bce117adbd24240ba026f705e782d185)), closes [nestjs/nest#5990](https://github.com/nestjs/nest/issues/5990) [#3349](https://github.com/biosimulations/biosimulations/issues/3349)
* **exceptions:** improve error handling ([083091e](https://github.com/biosimulations/biosimulations/commit/083091e8103340d6fd80534ec05cbc0d815e81cf))


### Features

* **api:** add caching to results and ontology endpoints ([bb2a991](https://github.com/biosimulations/biosimulations/commit/bb2a991c77ee32bb43ae9edb02938ca345700544))
* **api:** add caching to results endpoints ([ed54363](https://github.com/biosimulations/biosimulations/commit/ed5436370f4a9be2ed253d69417144d0db10e89b))
* **api:** add health check for job queue ([a05556e](https://github.com/biosimulations/biosimulations/commit/a05556ef33817253b32e52e8519806599e966723))
* **api:** add health module and endpoints ([da036f5](https://github.com/biosimulations/biosimulations/commit/da036f59d5f2479790a5d4827a2b830432ea35af))
* **api:** add various health checks and endpoints ([d30c7d2](https://github.com/biosimulations/biosimulations/commit/d30c7d25d62268f675cff1804b38309747c68e0b))
* **api:** setup results cache with REDIS ([3e40fde](https://github.com/biosimulations/biosimulations/commit/3e40fdebe38ca3e3765823bf7f98fd33895a7c9f))
* **dispatch-service:** add limit for retries of status for jobs ([f187b03](https://github.com/biosimulations/biosimulations/commit/f187b03024215b0cb2d53858085197b9d45736b9))
* **dispatch,platform:** added structured data for projects, simulation runs: ([7135ed0](https://github.com/biosimulations/biosimulations/commit/7135ed01f3107407906b81ad5b9ac94846f36f82))
* **dispatch,simulators:** added structured data tutorials ([5c23e8d](https://github.com/biosimulations/biosimulations/commit/5c23e8da1bcc2e12163f82aebf98c3a90c14bab8))
* **dispatch,simulators:** encoded FAQs into Schema.org ([936fb44](https://github.com/biosimulations/biosimulations/commit/936fb4499a2b132945b39fc2148428d26a268189))
* **exceptions:** dont process health check http exceptions ([b75747e](https://github.com/biosimulations/biosimulations/commit/b75747e6e04c511899037322e4829ae12433f3ab))
* **simulators-api,api:** added clearer payload too large messages ([98fb49d](https://github.com/biosimulations/biosimulations/commit/98fb49d651c01d621d175051cb030621b273034a))
* **simulators-api:** add health checks for simulators-api ([26d3b63](https://github.com/biosimulations/biosimulations/commit/26d3b63d37b1f271c0fc327b4a8e2c4981565651))
* **simulators:** added structured data for simulators as software applications ([7618db3](https://github.com/biosimulations/biosimulations/commit/7618db3ee539c31be2167a8d717cb4a5fac2c798))


### Reverts

* **simulators-api,api:** revert partially the changes in 98fb49d651c01d621d175051cb030621b273034a ([d88fcea](https://github.com/biosimulations/biosimulations/commit/d88fcead0eebe59f6394c9befca9e6ac7e132c83))

## [5.8.0](https://github.com/biosimulations/biosimulations/compare/v5.7.3...v5.8.0) (2021-10-20)


### Features

* **api:** enabling simulation run requests with latest version of a simulator ([2b7c2d9](https://github.com/biosimulations/biosimulations/commit/2b7c2d92ed10825bfce1a6c35a5a1b4908ceeebe))
* **config:** added endpoint for latest versions of simulators ([2a582ee](https://github.com/biosimulations/biosimulations/commit/2a582ee364a6308f42ae632b2487edd050964c87))
* **dispatch:** recorded simulator versions and digests for simulation runs ([f2abb39](https://github.com/biosimulations/biosimulations/commit/f2abb39790aa24c582c42521aca38f3dfbecaa56))
* **simulators:** added validation that version isn't reserved word 'latest' ([47843a8](https://github.com/biosimulations/biosimulations/commit/47843a8b687ff2ad001de286a07aa813735ca0c2))

## [5.7.3](https://github.com/biosimulations/biosimulations/compare/v5.7.2...v5.7.3) (2021-10-19)


### Bug Fixes

* **api:** correct field name to get values from dataservice ([53a6bbc](https://github.com/biosimulations/biosimulations/commit/53a6bbc5d2ccab9ea086f405656bf26d0cb2bacb)), closes [#3313](https://github.com/biosimulations/biosimulations/issues/3313)
* **api:** fix typo with checks ([caada0a](https://github.com/biosimulations/biosimulations/commit/caada0ab6e0f806508ee4d29966e5bcb0ca85109))

## [5.7.2](https://github.com/biosimulations/biosimulations/compare/v5.7.1...v5.7.2) (2021-10-19)


### Bug Fixes

* **api:** aligned parameter name in documentation ([0100a27](https://github.com/biosimulations/biosimulations/commit/0100a274725e22654ce4cc4534e5b0ff64a42de1))
* **simulators-api:** corrected put method; closes [#3305](https://github.com/biosimulations/biosimulations/issues/3305) ([57c34be](https://github.com/biosimulations/biosimulations/commit/57c34bef1f1370fdcd0a235020e4e7b5f20f5d54))

## [5.7.1](https://github.com/biosimulations/biosimulations/compare/v5.7.0...v5.7.1) (2021-10-18)


### Bug Fixes

* **account-api:** fixed route parameter names ([c7d8712](https://github.com/biosimulations/biosimulations/commit/c7d8712b0614d325bc1bf7fcc9145effcb6f2e60))
* **api:** fixed route parameter names ([1496f0f](https://github.com/biosimulations/biosimulations/commit/1496f0fc1c9ca0d63a620f73c216d9ab0752cce2))

## [5.7.0](https://github.com/biosimulations/biosimulations/compare/v5.6.2...v5.7.0) (2021-10-18)


### Bug Fixes

* **api:** fix docs and typing of open api definition ([#3307](https://github.com/biosimulations/biosimulations/issues/3307)) ([0640c6a](https://github.com/biosimulations/biosimulations/commit/0640c6aff0dd6a2d558d44be620d042b6e7ba49d)), closes [#3304](https://github.com/biosimulations/biosimulations/issues/3304)
* **api:** fix param name parsing for projectId ([af3c406](https://github.com/biosimulations/biosimulations/commit/af3c4069f352839a33696e3cf7c5dbb45d20c210))
* **combine-service:** added missing Swagger templates to Docker image ([f27b8c8](https://github.com/biosimulations/biosimulations/commit/f27b8c8831b11fffedd3355bc9668f78ad2e080c))


### Features

* **combine-service:** added health endpoint ([0d356d3](https://github.com/biosimulations/biosimulations/commit/0d356d347f11c3377f14ddbb66af8823198eaa30))
* **combine-service:** increased file upload limit, clarified error message ([5dfa25c](https://github.com/biosimulations/biosimulations/commit/5dfa25cca9fe73e2ef4b1af881f646bf9b022d5e))

## [5.6.2](https://github.com/biosimulations/biosimulations/compare/v5.6.1...v5.6.2) (2021-10-18)


### Bug Fixes

* **api:** fix module import ([d47c376](https://github.com/biosimulations/biosimulations/commit/d47c376915713fbc2ae97b296318061374b4bc10))
* **dispatch:** corrected file types for validate OMEX metadata form ([98794cd](https://github.com/biosimulations/biosimulations/commit/98794cd1397cb1ee921fad02b09d83ca82e7bd3f))

## [5.6.1](https://github.com/biosimulations/biosimulations/compare/v5.6.0...v5.6.1) (2021-10-18)


### Bug Fixes

* **api:** fix permissions for endpoints ([f00f6d1](https://github.com/biosimulations/biosimulations/commit/f00f6d116cc64525656f4b36030ee0109f9ff3b0)), closes [#3242](https://github.com/biosimulations/biosimulations/issues/3242)
* update client ids for api docs ([1ec36bb](https://github.com/biosimulations/biosimulations/commit/1ec36bb377ab939062da27455d199dbc3e4ada25))
## [5.6.0](https://github.com/biosimulations/biosimulations/compare/v5.5.0...v5.6.0) (2021-10-17)


### Bug Fixes

* **api:** update api to use updated hsds client ([e7832aa](https://github.com/biosimulations/biosimulations/commit/e7832aaf2c6d34c78903227dfc7e7db254b5b73d))
* **dispatch:** add flag to skip downloading test results to simulators service ([cbd63cc](https://github.com/biosimulations/biosimulations/commit/cbd63ccc6ecf757dec0da6cf000c200f2e8ebc4c)), closes [#3197](https://github.com/biosimulations/biosimulations/issues/3197)
* **dispatch:** fix alg list to empty list to prevent crashing ([338be99](https://github.com/biosimulations/biosimulations/commit/338be997320e053b054137a44546a82a21c13e34))
* make changes to update mongoose ([4653155](https://github.com/biosimulations/biosimulations/commit/46531556323e5d84cc96afbd87f4a58aa335be23))
* **platform,dispatch:** corrected link to simulation results in files tab ([00aa363](https://github.com/biosimulations/biosimulations/commit/00aa363a059fc9b0f0f95d316c1313efadb890b7))
* **platform:** redirect to 404 for non-existent projects; closes [#3234](https://github.com/biosimulations/biosimulations/issues/3234) ([173439a](https://github.com/biosimulations/biosimulations/commit/173439a2c322d933dc2ad4f5f1b3bcbeee666e80))
* **simulators:** fix json-ld metadata on index.html ([95b3d98](https://github.com/biosimulations/biosimulations/commit/95b3d983c1be2ddf633570214f8d35a920d98f1e))


### Features

* **api:** add custom styling to swagger ui ([90830c0](https://github.com/biosimulations/biosimulations/commit/90830c02f53df0fc9f86c7211b2fa3359a1b8807))
* **hsds:** update client library ([5746832](https://github.com/biosimulations/biosimulations/commit/574683246e7abb8a338d35d2a7a98ba25ce96d41))
* **simulators:** added repository digest to image model; closes [#3194](https://github.com/biosimulations/biosimulations/issues/3194) ([1293410](https://github.com/biosimulations/biosimulations/commit/129341077889888b41397cf471389d8921c8c2f8))
* **simulators:** expanded full text search; closes [#3209](https://github.com/biosimulations/biosimulations/issues/3209) ([c877189](https://github.com/biosimulations/biosimulations/commit/c877189fc5ece4943f3577602ff770660cdf01c0))


### Reverts

* **deps:** revert update dependency eslint to v8 ([0fdb3d8](https://github.com/biosimulations/biosimulations/commit/0fdb3d81a59710e4667974ecf0e42c8ec65ffd34))

## [5.5.0](https://github.com/biosimulations/biosimulations/compare/v5.4.0...v5.5.0) (2021-10-09)


### Bug Fixes

* **simulators-api:** add biosimulations.org to cors ([c2eea89](https://github.com/biosimulations/biosimulations/commit/c2eea89c81f266a00235bd1338102bbee7274dae))


### Features

* **api:** add case-insenstive unique index for project ids ([5f96f91](https://github.com/biosimulations/biosimulations/commit/5f96f91d83ae8835a313d6273eafab5e544f4e77)), closes [#3160](https://github.com/biosimulations/biosimulations/issues/3160)
* **api:** add controller level validation for project ids ([cdba9ce](https://github.com/biosimulations/biosimulations/commit/cdba9ce007fa92515f64147527824ca3df225808))
* **dispatch:** added check that simulation run was successful ([b4ade32](https://github.com/biosimulations/biosimulations/commit/b4ade32f2c3443d182f92ced5d12bf80341a260a))
* **platform:** added validation for project ids; closes [#3183](https://github.com/biosimulations/biosimulations/issues/3183) ([01b6178](https://github.com/biosimulations/biosimulations/commit/01b61788565e3c42c1b11acf1e2f02b349321cfc))

## [5.4.0](https://github.com/biosimulations/biosimulations/compare/v5.3.0...v5.4.0) (2021-10-08)


### Bug Fixes

* **combine-service:** update combine-service client to latest api changes ([74a70d8](https://github.com/biosimulations/biosimulations/commit/74a70d8fa73086277318b363dc58d1ebfa7d1970))
* **dispatch,platform:** fixed visualization rendering when groups of datasets are selected ([f205ab4](https://github.com/biosimulations/biosimulations/commit/f205ab457eab7712133d859e6ae4cfeb8dd16d63))
* **dispatch:** handle cases when metadata is empty without throwing error ([57a16c5](https://github.com/biosimulations/biosimulations/commit/57a16c54e496929a482fc69edc9bf1b2ca165e9d))


### Features

* **combine-service:** added endpoints for validation ([f602e65](https://github.com/biosimulations/biosimulations/commit/f602e651107f7efd195bb05ef41650bfe64009cd))
* **combine-service:** updated to biosimulators-utils 0.1.130 ([e46d13e](https://github.com/biosimulations/biosimulations/commit/e46d13ecc28e0e86ff2377068f8e11872276d327))
* **dispatch:** add some error handling ([6cc7f62](https://github.com/biosimulations/biosimulations/commit/6cc7f625451306c04e53884ab675d33ffd1fd5b8)), closes [#3088](https://github.com/biosimulations/biosimulations/issues/3088)
* **dispatch:** added forms for validating models, simulations and metadata ([d7991ed](https://github.com/biosimulations/biosimulations/commit/d7991ed6ce4d5cc14e35fbd70143ebb1d21705ca))
* **dispatch:** added options to project validation ([eebd16f](https://github.com/biosimulations/biosimulations/commit/eebd16fb26a66515952708e1f8e11638f234e93e))


### Reverts

* **combine-service:** reverted URL for COMBINE API ([d1dd8b5](https://github.com/biosimulations/biosimulations/commit/d1dd8b595ae1777c4e73b1e90d83f6a33603730c))

## [5.3.0](https://github.com/biosimulations/biosimulations/compare/v5.2.0...v5.3.0) (2021-10-07)


### Bug Fixes

* added 'master' attribute to file object ([5f4f722](https://github.com/biosimulations/biosimulations/commit/5f4f722d848e1e0a37802119b860b945677417a6))
* **api,mail-service:** update api client to return observable ([d26ce89](https://github.com/biosimulations/biosimulations/commit/d26ce8906de9275f40790ea4afbe6a97e404dd0d)), closes [#3102](https://github.com/biosimulations/biosimulations/issues/3102)
* **api:** add permissions to get all specs ([dbce421](https://github.com/biosimulations/biosimulations/commit/dbce4210625bc6ed6bb64eb7e3a170e385c88e55)), closes [#3136](https://github.com/biosimulations/biosimulations/issues/3136)
* **deps:** update dependency rxjs to v7.3.1 ([8cb2d32](https://github.com/biosimulations/biosimulations/commit/8cb2d3209fbcdc58771c4d5517aa5095376e87b5))
* **deps:** update nest ([f7a97e6](https://github.com/biosimulations/biosimulations/commit/f7a97e60127d9e66103dd331ce51038431398433))
* **dispatch:** fixed lint issue ([660bcb8](https://github.com/biosimulations/biosimulations/commit/660bcb8dcb6957bdcfc0d7469fcafb21f29f84c9))
* fixed lint issue ([e82e6fc](https://github.com/biosimulations/biosimulations/commit/e82e6fc06df38cc330dd7539e0cc7e18c45caab1))
* **platform,dispatch:** fixed plotly tests ([f2d8353](https://github.com/biosimulations/biosimulations/commit/f2d83535195d6d1716f2c3e9c3dc5d74aa310852))
* **simulators-api:** allow cors for biosimulatiors.dev ([1d93452](https://github.com/biosimulations/biosimulations/commit/1d93452b12174b444c5d53417521660dde72be56))
* **ui:** fixed display of errors with Vega visualizations ([1270955](https://github.com/biosimulations/biosimulations/commit/12709555792627db191d6794477f72ca6d81c7c4))
* **ui:** fixed Vega export for 1-d heatmaps ([cd5713d](https://github.com/biosimulations/biosimulations/commit/cd5713d73f14cddac1febaeb977ebb77a5fa0dff))


### Features

* **api:** create project endpoints ([d1b9fe7](https://github.com/biosimulations/biosimulations/commit/d1b9fe73c719358ace007c9d67a43c4a1d1c6810)), closes [#3067](https://github.com/biosimulations/biosimulations/issues/3067)
* **combine-service:** added options for validation of COMBINE archives ([42febbe](https://github.com/biosimulations/biosimulations/commit/42febbedc289d76f41bd1655014e1c9a172ff643))
* **combine-service:** added options to control COMBINE archive validation ([b4c0c12](https://github.com/biosimulations/biosimulations/commit/b4c0c123ea6ac776cc912e4eed964f773477d403))
* **combine-service:** added timeout for simulation execution ([8eb8deb](https://github.com/biosimulations/biosimulations/commit/8eb8debcab97b215260b9e6646751850b1a47593))
* **combine-service:** update combine-api client ([77c2f6d](https://github.com/biosimulations/biosimulations/commit/77c2f6df36cc59be74d4271938c34b5f074608a1))
* **datamodel:** add project datamodel ([57cf45c](https://github.com/biosimulations/biosimulations/commit/57cf45ceda0ba7d32f909e5d07fb7200e8dbbdee))
* **dispatch,platform,simulators:** improve recognition of Vega files by media type ([5f0051b](https://github.com/biosimulations/biosimulations/commit/5f0051bbfbd4e2ab6cd51d98a4352487ea94d640))
* **dispatch:** added options for validating COMBINE archives ([7d0a815](https://github.com/biosimulations/biosimulations/commit/7d0a8159798132fd65d93d7c91881a834b389c83))
* **platform:** added export of visualizations to Vega and COMBINE archives ([fce6731](https://github.com/biosimulations/biosimulations/commit/fce67312fa9432ed33ee0714a30d81d1496de111))
* **platform:** added heatmap, line plots ([84898db](https://github.com/biosimulations/biosimulations/commit/84898db3682a5b9321132a8a3058a6663c56113a))
* **platform:** added histogram visualization ([7a6abfa](https://github.com/biosimulations/biosimulations/commit/7a6abfafe6fad3dcebccff2977fec994839cdcbc))
* **platform:** added SED-ML visualizations ([2ac40b3](https://github.com/biosimulations/biosimulations/commit/2ac40b37cf40c713fd798c92fa761aa080586077))
* **platform:** added simulation types and algorithms to simulaton overview ([4c768e5](https://github.com/biosimulations/biosimulations/commit/4c768e594caee33d7cee79f157eaa70ba6703a3c))
* **platform:** added vega export for heatmaps and line plots ([a8c9ad4](https://github.com/biosimulations/biosimulations/commit/a8c9ad418407d7ec7d7f42e21ddf6a0d0dd5d966))
* **platform:** front end displays projects from api ([9ecfa80](https://github.com/biosimulations/biosimulations/commit/9ecfa80ce7f70a45762563b88c82c0bf0e0cf3a0)), closes [#3149](https://github.com/biosimulations/biosimulations/issues/3149)
* **ui:** added ability to attach hyperlinks to menu items ([7ca3d10](https://github.com/biosimulations/biosimulations/commit/7ca3d10209fa04d77cbc220664f0d1870c542c12))


### Reverts

* **deps:** revert 235c9db3e9649cdb8b42e6575517aa651f9e1c2d ([05cb6f3](https://github.com/biosimulations/biosimulations/commit/05cb6f39142d576e4d9c866a6ee2177589aadbd3))

## [5.2.0](https://github.com/biosimulations/biosimulations/compare/v5.1.1...v5.2.0) (2021-10-04)


### Bug Fixes

* **api:** add authentication to post metadata ([999e2b9](https://github.com/biosimulations/biosimulations/commit/999e2b9e11ddaa26e001d89ccbbeb26541cc3f5a)), closes [#2865](https://github.com/biosimulations/biosimulations/issues/2865)
* **deps:** update dependency @sendgrid/mail to v7.4.7 ([bae8e96](https://github.com/biosimulations/biosimulations/commit/bae8e968c0062da5f477ca514e3ed3c08e6686e3))
* **dispatch-service:** fix processing of environment variables ([cb0aa04](https://github.com/biosimulations/biosimulations/commit/cb0aa04cbdc7bd8f2a3e7cf818c541470d7c4519))
* **dispatch-service:** process metadata with other processing ([7f8f44a](https://github.com/biosimulations/biosimulations/commit/7f8f44aecc741dcd746e0de71489bcd5a9c68319)), closes [#3046](https://github.com/biosimulations/biosimulations/issues/3046)


### Features

* **api:** add endpoints to get particular specifications for simulation runs ([87eede1](https://github.com/biosimulations/biosimulations/commit/87eede1c0283bdc3a4a5b6614dba06354a96ddf6))
* **api:** add file object ([39f25f3](https://github.com/biosimulations/biosimulations/commit/39f25f3668d7e249ae1b9834b6ef66cb58350208)), closes [#2914](https://github.com/biosimulations/biosimulations/issues/2914)
* **api:** create specifications object and endpoints ([aca4786](https://github.com/biosimulations/biosimulations/commit/aca4786282526813eab3280c87d495617c7a2ef8))
* **dispatch-service:** process files and sedml specs ([fb37624](https://github.com/biosimulations/biosimulations/commit/fb376243c33043d4945a2849d0df50a54332fcda))
* **dispatch-service:** send sedml specifications to the api ([7d6a8c7](https://github.com/biosimulations/biosimulations/commit/7d6a8c7baec360f985459534f9fd5d67e4342260))

## [5.1.1](https://github.com/biosimulations/biosimulations/compare/v5.1.0...v5.1.1) (2021-10-01)


### Bug Fixes

* **platform:** corrected handling of software license keys for simulation ([84bef20](https://github.com/biosimulations/biosimulations/commit/84bef20d9c035168d7d74fcffdeafae97091af8f))
* **platform:** corrected handling of software license keys for simulation ([aba13fe](https://github.com/biosimulations/biosimulations/commit/aba13fe9c30c2355d6d63a9a3546fa78062b82d6))

## [5.1.0](https://github.com/biosimulations/biosimulations/compare/v5.0.0...v5.1.0) (2021-09-29)


### Bug Fixes

* **deps:** update dependency auth0 to v2.36.2 ([#3076](https://github.com/biosimulations/biosimulations/issues/3076)) [skip ci] ([a696fbc](https://github.com/biosimulations/biosimulations/commit/a696fbcc9f495fef5b20fba7235a8fb7d590b4d9))
* **dispatch-api:** fix cors for biosimulators ([22439d2](https://github.com/biosimulations/biosimulations/commit/22439d2d6d2c16ce71dc652e5353a9f902a29c6a))


### Features

* **dispatch:** require configuration of academic use for commercial solvers ([6c5307c](https://github.com/biosimulations/biosimulations/commit/6c5307c31dad3f59b2e7c01b07ea1b83218a7ff0))

## [5.0.0](https://github.com/biosimulations/biosimulations/compare/v4.6.0...v5.0.0) (2021-09-28)


### Bug Fixes

* **dispatch-service:** fix type error when processing metadata ([a646c98](https://github.com/biosimulations/biosimulations/commit/a646c98ef2433b5b20b754fd637f326513aa57e1))
* make datamodel consistent for license ([4b95e4d](https://github.com/biosimulations/biosimulations/commit/4b95e4d89af83b69334602f53ffacaa0744e5aff)), closes [#3050](https://github.com/biosimulations/biosimulations/issues/3050)
* **dispatch-api:** remove extra slash for metadata uris ([b627e74](https://github.com/biosimulations/biosimulations/commit/b627e7491995ae1c1e70feea93b1e7f4cc53902a)), closes [#3052](https://github.com/biosimulations/biosimulations/issues/3052)
* ensure external url is used for combine api ([2d98aba](https://github.com/biosimulations/biosimulations/commit/2d98aba10be03163e9b43aa67da69a95910eb763))
* **dispatch:** corrected when metadata about simulation projects is retrieved ([1682d83](https://github.com/biosimulations/biosimulations/commit/1682d83d6870e4bcfa1d4b895d3a88d2cee60285))


### Code Refactoring

* consolidate backend apis ([b27bd0e](https://github.com/biosimulations/biosimulations/commit/b27bd0e260336df3553b1b3a7e3447c0e26ac716)), closes [#2724](https://github.com/biosimulations/biosimulations/issues/2724)


### Features

* **dispatch-api:** ensure only public models are shown for platform ([#3045](https://github.com/biosimulations/biosimulations/issues/3045)) ([5619c03](https://github.com/biosimulations/biosimulations/commit/5619c03d7fc088d0b3be33136935c80e7cb9c862)), closes [#3044](https://github.com/biosimulations/biosimulations/issues/3044)
* **dispatch-api:** extract files to s3 and replace combine archive file extraction endpoint ([56f8413](https://github.com/biosimulations/biosimulations/commit/56f84133b193c4d54f77c33ba2c01105df6162e3)), closes [#2945](https://github.com/biosimulations/biosimulations/issues/2945)
* **dispatch-api:** upload omex files to s3 from url ([4d8f780](https://github.com/biosimulations/biosimulations/commit/4d8f78058d3ea7050d913ad651c907c0a631a3f4))


### Reverts

* **dispatch-api:** revert permissions change in 3175c6378160f34e8389b6e501ea2534eb9d4c12 ([5ab7d08](https://github.com/biosimulations/biosimulations/commit/5ab7d08cc922778436646dca0331aad7bafef0d3))


### BREAKING CHANGES

* The ontology, dispatch, and platform apis are consolidated into one main backend
api for biosimulations. There is a seperate api for biosimulators. The combine-service also provides a rest api that is mostly intended for internal use.

## [4.6.0](https://github.com/biosimulations/biosimulations/compare/v4.5.0...v4.6.0) (2021-09-27)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.995.0 ([19b509a](https://github.com/biosimulations/biosimulations/commit/19b509a7100e59dbcdf3c2674262ec2bf4333f77))
* **dispatch:** improved handling of undefined simulations in simulations browse view; closes [#2999](https://github.com/biosimulations/biosimulations/issues/2999) ([bf47994](https://github.com/biosimulations/biosimulations/commit/bf4799404320d07437797804973e933d7f147c7e))
* **platform:** handle cases when metadata is missing ([f738a74](https://github.com/biosimulations/biosimulations/commit/f738a741caba0da783e522c8d74e86c04b9e4aa8))


### Features

* **simulators:** improved specification of simulator CLIs; closes [#3015](https://github.com/biosimulations/biosimulations/issues/3015) ([c396bfa](https://github.com/biosimulations/biosimulations/commit/c396bfa79a7f9abbc0a3e6b374f96d726fe5eaa7))
* create alternate vega view component ([453a109](https://github.com/biosimulations/biosimulations/commit/453a1094677a18ccc2019bd953d618e4e81c6e2f))
* **dispatch:** added license confirmation to publish run form ([1889e1e](https://github.com/biosimulations/biosimulations/commit/1889e1ec4f50579c51838494f7d693fdb5a87415))
* **dispatch,platform:** updated terms about granting BioSimulations to distribute projects ([6ee2bb1](https://github.com/biosimulations/biosimulations/commit/6ee2bb1779e3f02df9dee82bde75e11fa5dc5f9d))
* **platform:** add code to get metadata and specs of visualiazations ([0cdc472](https://github.com/biosimulations/biosimulations/commit/0cdc4727b201858adec4cb40565bb643c01c2905))
* **platform:** add support for showing vega figures ([f011480](https://github.com/biosimulations/biosimulations/commit/f011480096ee29f42e8aab833dbaec89a1bc1596))
* **ui:** allow for conditional loading of tabs ([1b4779f](https://github.com/biosimulations/biosimulations/commit/1b4779f177748922805e8aaad196d0ef412a5239))
* updated biosimulators-utils, biosimulators-bionetgen ([ce75ea6](https://github.com/biosimulations/biosimulations/commit/ce75ea6a46e880b1c021a336992d7acc311f200d))


### Performance Improvements

* **platform:** fix tests ([a8f6d0e](https://github.com/biosimulations/biosimulations/commit/a8f6d0ed968a9ccccd45249db4b71313ef7a5f5b))

## [4.5.0](https://github.com/biosimulations/biosimulations/compare/v4.4.2...v4.5.0) (2021-09-26)


### Features

* **simulators:** improved simulator usage examples; closes [#3016](https://github.com/biosimulations/biosimulations/issues/3016) ([a56dea2](https://github.com/biosimulations/biosimulations/commit/a56dea24e33e33b7db215a0cc1f62f4bc94dac7d))

## [4.4.2](https://github.com/biosimulations/biosimulations/compare/v4.4.1...v4.4.2) (2021-09-23)


### Bug Fixes

* **simulators-api:** fixed sorting of simulator versions with > 4 points ([020073a](https://github.com/biosimulations/biosimulations/commit/020073ab0a9e8eb0344d7ab4ffee287a323caaed)), closes [#3008](https://github.com/biosimulations/biosimulations/issues/3008)


### Reverts

* **simulators-api:** revert [#3009](https://github.com/biosimulations/biosimulations/issues/3009) to prevent container crashing ([c7ff2cc](https://github.com/biosimulations/biosimulations/commit/c7ff2cc018e28dad55f771c6e778193082e9c2ff)), closes [#3008](https://github.com/biosimulations/biosimulations/issues/3008)

## [4.4.1](https://github.com/biosimulations/biosimulations/compare/v4.4.0...v4.4.1) (2021-09-22)


### Bug Fixes

* fixed sorting of simulator versions with > 4 points ([e0b60ce](https://github.com/biosimulations/biosimulations/commit/e0b60ce7e857a949d73f4daed4cf51e9e2b0ca91))

## [4.4.0](https://github.com/biosimulations/biosimulations/compare/v4.3.0...v4.4.0) (2021-09-22)


### Features

* **ontology:** updated to KiSAO 2.29 ([31e7d7b](https://github.com/biosimulations/biosimulations/commit/31e7d7ba7115c9f381b636113221edca9010397b))
* **platform:** improve styling of platform browse ([28e8e1d](https://github.com/biosimulations/biosimulations/commit/28e8e1d6858363325ae29ae3d71e1dea2a1b19c9))


### Reverts

* remove ui commit scope [skip ci] ([5051250](https://github.com/biosimulations/biosimulations/commit/50512504bd98cf55f2111ad0bf074ae5837260ea))

## [4.3.0](https://github.com/biosimulations/biosimulations/compare/v4.2.0...v4.3.0) (2021-09-15)


### Bug Fixes

* **deps:** update dependency @stoplight/json-ref-resolver to v3.1.3 ([#2986](https://github.com/biosimulations/biosimulations/issues/2986)) ([fe9c5f3](https://github.com/biosimulations/biosimulations/commit/fe9c5f372c69c2f8a11d0d401a312ecdb7bf3338))
* **deps:** update dependency bull to v3.29.2 ([#2987](https://github.com/biosimulations/biosimulations/issues/2987)) ([98cbebf](https://github.com/biosimulations/biosimulations/commit/98cbebf7752e6fd6b034e0a459b6abc9de106247))
* **dispatch:** proceed if metadata is missing ([#2998](https://github.com/biosimulations/biosimulations/issues/2998)) ([f09c633](https://github.com/biosimulations/biosimulations/commit/f09c633a4b69b096c5bd07a1377c451ea3ae3aa1)), closes [#2994](https://github.com/biosimulations/biosimulations/issues/2994)


### Features

* **simulators:** expanded specs for simulators ([32b100b](https://github.com/biosimulations/biosimulations/commit/32b100be3f5856bc8131427155031dc5abe1013a))
* expanded simulator specs ([f281cd3](https://github.com/biosimulations/biosimulations/commit/f281cd31b7c3c7c08dfc47162944dcfdbb7c4761))

## [4.2.0](https://github.com/biosimulations/biosimulations/compare/v4.1.0...v4.2.0) (2021-09-11)


### Bug Fixes

* **deps:** update dependency axios to v0.21.4 ([#2952](https://github.com/biosimulations/biosimulations/issues/2952)) ([c24c0e2](https://github.com/biosimulations/biosimulations/commit/c24c0e2eb96d6ff5c1a77f9408dee9c80e02c0f1))
* **deps:** update dependency ssh2 to v1.4.0 ([#2953](https://github.com/biosimulations/biosimulations/issues/2953)) ([a3afb2b](https://github.com/biosimulations/biosimulations/commit/a3afb2b7c862b7cf4bb86a451dce7380a64afa39))
* specify global setTimeout instead of window ([432b90a](https://github.com/biosimulations/biosimulations/commit/432b90a00f34c604bc5f2fd7038a925c52ea4fec))
* **dispatch-service:** restore check for empty env variables ([3f714ba](https://github.com/biosimulations/biosimulations/commit/3f714baf82ab04673fb33a28cc9f7daa9899c39b))
* restore mkdocs file location ([41a269e](https://github.com/biosimulations/biosimulations/commit/41a269e37f9c3e704baa614b4ab8de30d8ed6546))


### Features

* **account-api:** replace typegoose with mongoose ([69911a3](https://github.com/biosimulations/biosimulations/commit/69911a35ab0ca6c324ba53a097c814d53730d491))
* **dispatch-api:** handle errors and timeouts on uploa> ([dad35ea](https://github.com/biosimulations/biosimulations/commit/dad35ea5718abac450c682970e332ee4292890f3)), closes [#2860](https://github.com/biosimulations/biosimulations/issues/2860)
* **dispatch-service:** added passing software licenses from deployment secrets to Singularity run ([cc19999](https://github.com/biosimulations/biosimulations/commit/cc199990c6255b6b70bc436a9d16051602c4d0c5))
* **storage:** add simulation storage service and timeout for s3 uploads ([0c24173](https://github.com/biosimulations/biosimulations/commit/0c241737989289f79b13ebc8efa265cc7c6fc91f))

## [4.1.0](https://github.com/biosimulations/biosimulations/compare/v4.0.1...v4.1.0) (2021-09-06)


### Features

* **simulators:** added attribute to track installation instructions for Python APIs ([cb2b415](https://github.com/biosimulations/biosimulations/commit/cb2b415b4513170bbb09140e6cd9bff4b970d3ed))

## [4.0.1](https://github.com/biosimulations/biosimulations/compare/v4.0.0...v4.0.1) (2021-09-06)


### Bug Fixes

* **dispatch:** simulation results URLs for data visualizations ([9b7b879](https://github.com/biosimulations/biosimulations/commit/9b7b8795ac2b524b9089fad5a2c7916e3d1214f4))


### Reverts

* 39a60b17d640b62639f6594024f4ba4c66baedc5 ([f804cce](https://github.com/biosimulations/biosimulations/commit/f804cce9e3b3787a11b2989743e86407a4c014dd)), closes [#2959](https://github.com/biosimulations/biosimulations/issues/2959)

## [4.0.0](https://github.com/biosimulations/biosimulations/compare/v3.20.0...v4.0.0) (2021-09-04)


### Bug Fixes

* **dispatch:** properly encode uri to allow for fetching results ([dcbf044](https://github.com/biosimulations/biosimulations/commit/dcbf04433f7e105a01f501f3aa7172c82807ea41))


### Features

* update example simulation runs ([395f513](https://github.com/biosimulations/biosimulations/commit/395f513657c662d2b26b3d3b0de95cdd860ea326)), closes [#2951](https://github.com/biosimulations/biosimulations/issues/2951)


### BREAKING CHANGES

* simulation runs sumbitted prior to the update will not display on the dispatch app

## [3.20.0](https://github.com/biosimulations/biosimulations/compare/v3.19.0...v3.20.0) (2021-09-04)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.983.0 ([#2947](https://github.com/biosimulations/biosimulations/issues/2947)) ([39a60b1](https://github.com/biosimulations/biosimulations/commit/39a60b17d640b62639f6594024f4ba4c66baedc5))
* **dispatch-service:** correct the determination of the environment ([ce46d3b](https://github.com/biosimulations/biosimulations/commit/ce46d3bddb05195dd29408d05c00de53186336ae))
* fix default environment to dev ([c970ccd](https://github.com/biosimulations/biosimulations/commit/c970ccde7a533bb0db3f0ec8334308d3a1ee237d))
* new endpoint implementation ([ed42b6b](https://github.com/biosimulations/biosimulations/commit/ed42b6b27fbba4a97b8a931d6d9771b1563eddb9)), closes [#2943](https://github.com/biosimulations/biosimulations/issues/2943) [#2861](https://github.com/biosimulations/biosimulations/issues/2861) [#2859](https://github.com/biosimulations/biosimulations/issues/2859)


### Features

* **dispatch,dispatch-api:** move thumbnail processing to backend ([4495d6d](https://github.com/biosimulations/biosimulations/commit/4495d6d70e8fcb168fd5b4a38f70850171908d7b))
* **platform:** add page to view projects on platform ([e568d0c](https://github.com/biosimulations/biosimulations/commit/e568d0c16a5e67198e921d8705e45f137c680df2))


### Reverts

* revert commit 5dad745d1df0ffc3fb2fba8fc3b99b21b69b0521 ([f8cdd5b](https://github.com/biosimulations/biosimulations/commit/f8cdd5b338ec3fd7f8b7faf607ce893a9d343075))

## [3.19.0](https://github.com/biosimulations/biosimulations/compare/v3.18.0...v3.19.0) (2021-09-02)


### Features

* **combine-service:** updated to biosimulators-utils 0.1.115, biosimulators-amici 0.1.18 ([9ad2945](https://github.com/biosimulations/biosimulations/commit/9ad29450a51b8ff181a00fe57c70b660dc917a60))

## [3.18.0](https://github.com/biosimulations/biosimulations/compare/v3.17.0...v3.18.0) (2021-09-01)


### Bug Fixes

* **deps:** update dependency form-data to v4 ([#2925](https://github.com/biosimulations/biosimulations/issues/2925)) ([36a79ba](https://github.com/biosimulations/biosimulations/commit/36a79baccabb13e181e9ca5ee0cd2d0bff629697))
* **deps:** update dependency jwks-rsa to v2 ([#2928](https://github.com/biosimulations/biosimulations/issues/2928)) ([f4a3f10](https://github.com/biosimulations/biosimulations/commit/f4a3f107f5e43831e4cf72ffc63fcaf67a0026e3))
* **deps:** update dependency ssh2 to v1 ([#2930](https://github.com/biosimulations/biosimulations/issues/2930)) ([11111c7](https://github.com/biosimulations/biosimulations/commit/11111c76a5c943741397d3110189ac0d5ee53a86))


### Features

* **combine-service:** fixed error handling for run sim, simplified run sim options ([5e63d49](https://github.com/biosimulations/biosimulations/commit/5e63d49eb5a1ad5c27ac09dc970093f04ff79980))
* **dispatch:** added support for new SBO modeling framework terms ([80ee759](https://github.com/biosimulations/biosimulations/commit/80ee759d6be92545b01f999c1a7c0630fa43f43d))

## [3.17.0](https://github.com/biosimulations/Biosimulations/compare/v3.16.0...v3.17.0) (2021-09-01)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.980.0 ([#2906](https://github.com/biosimulations/Biosimulations/issues/2906)) ([163191d](https://github.com/biosimulations/Biosimulations/commit/163191d5cc1e24dbeb0440681761e175b633a759))


### Features

* add shared config file support ([976e578](https://github.com/biosimulations/Biosimulations/commit/976e57846a8c43fa10f8be4e70a8a1989bde683c))
* **dispatch:** call the metadata endpoint to get simulation metadata ([ae1054f](https://github.com/biosimulations/Biosimulations/commit/ae1054f6f101b170cb1408d13ffdcbb39f0b25a1)), closes [#2866](https://github.com/biosimulations/Biosimulations/issues/2866)

## [3.16.0](https://github.com/biosimulations/Biosimulations/compare/v3.15.0...v3.16.0) (2021-08-31)


### Bug Fixes

* **deps:** update dependency class-validator to v0.13.1 ([#2894](https://github.com/biosimulations/Biosimulations/issues/2894)) ([3676e59](https://github.com/biosimulations/Biosimulations/commit/3676e59d26a91d81bd7b12cd0287d619fa8e89ab))
* **deps:** update dependency nats to v2.2.0 ([#2886](https://github.com/biosimulations/Biosimulations/issues/2886)) ([7097edc](https://github.com/biosimulations/Biosimulations/commit/7097edc10eef06edc70ede51b1c2c4dc9eec5810))
* **deps:** update dependency rxjs to v7.3.0 ([#2895](https://github.com/biosimulations/Biosimulations/issues/2895)) ([c01604b](https://github.com/biosimulations/Biosimulations/commit/c01604bfcf4df5516158b442435f63002ef9720c))
* **deps:** update dependency stackdriver-errors-js to v0.10.0 ([43451aa](https://github.com/biosimulations/Biosimulations/commit/43451aa19accdc7c670a17cd6fec1aa660934234))
* **deps:** update dependency tslib to v2.3.1 ([#2888](https://github.com/biosimulations/Biosimulations/issues/2888)) ([fc78756](https://github.com/biosimulations/Biosimulations/commit/fc787565990a70489c6ee0a6b9450af7e92eb118))


### Features

* **dispatch:** added dry run option to example simulation submission ([1487880](https://github.com/biosimulations/Biosimulations/commit/1487880a87f82f48e29dafcfe9741bf9ff862cb7))
* **dispatch:** added example simulation run for RBApy ([277772e](https://github.com/biosimulations/Biosimulations/commit/277772e9d062f4a82fab455efbc2c67d088770ad))
* **dispatch-api:** set uris for metadata elements ([96e94fe](https://github.com/biosimulations/Biosimulations/commit/96e94fef4ea1a25ad95cba100c9b635618d29e7e))
* **simulators:** added ability to capture Python APIs in simulator specs ([5ed44cb](https://github.com/biosimulations/Biosimulations/commit/5ed44cb904349e4f91eb6b8cf62f264eb6eeebdd))

## [3.15.0](https://github.com/biosimulations/Biosimulations/compare/v3.14.0...v3.15.0) (2021-08-29)


### Bug Fixes

* **deps:** update dependency aws-sdk to v2.978.0 ([#2883](https://github.com/biosimulations/Biosimulations/issues/2883)) ([ca487aa](https://github.com/biosimulations/Biosimulations/commit/ca487aab71671f9a4a71668e253800fdc7e98708))
* **deps:** update dependency bull to v3.29.1 ([#2884](https://github.com/biosimulations/Biosimulations/issues/2884)) ([e6589e7](https://github.com/biosimulations/Biosimulations/commit/e6589e7d70baa6c9d6797d4fde95687d4798c19b))
* **deps:** update nest ([#2835](https://github.com/biosimulations/Biosimulations/issues/2835)) ([0e65b60](https://github.com/biosimulations/Biosimulations/commit/0e65b6084e75aae31ff8d089e6321f581fd8742d))


### Features

* **combine-service:** updated to Biosimulators-utils with support for RBA models ([610225b](https://github.com/biosimulations/Biosimulations/commit/610225b46f8bd9ed26c1ca632f01c051d5765dc8))
* updated biosimulators documentation links to docs.biosimulatos.org ([bfa49bb](https://github.com/biosimulations/Biosimulations/commit/bfa49bb7530cb656689eb8632b365110fb6b5aca))
* updated SBO for term for RBA ([fc64418](https://github.com/biosimulations/Biosimulations/commit/fc64418993cb06a162884461646b586801b7f37e))

## [3.14.0](https://github.com/biosimulations/Biosimulations/compare/v3.13.0...v3.14.0) (2021-08-27)


### Bug Fixes

* **combine-service:** dont change field "abstract" to "_abstract" ([591b0db](https://github.com/biosimulations/Biosimulations/commit/591b0db5bef0c4a260fb6f0583ed82b0eaf77481))
* **combine-service:** fix api client implementation ([d313f1c](https://github.com/biosimulations/Biosimulations/commit/d313f1cae8f64d1e93a0741bab2545d9b14d1109))
* **dispatch-service:** dont log error if job is not yet present ([b80db71](https://github.com/biosimulations/Biosimulations/commit/b80db71ed633cbf9c6cb2103a513d7dd1fc397cd))
* **dispatch-service:** fix url for posting metadata ([a5926e7](https://github.com/biosimulations/Biosimulations/commit/a5926e70d24985d0d50adfc65acdd4d6bbb8ca7d))
* **platform:** correct url for metadata ([c2a7a63](https://github.com/biosimulations/Biosimulations/commit/c2a7a63f57072598f9890905697edcfb461742f4))
* **platform:** fix unterminated string literal ([fe0343d](https://github.com/biosimulations/Biosimulations/commit/fe0343dc6fe808023c822b3ffeb004c094be96a2))


### Features

* **combine-service:** isolating simulation execution into separate processes ([9d23a5d](https://github.com/biosimulations/Biosimulations/commit/9d23a5d12571708ecf7b44bc83675cb5b5802a98))
* **combine-service:** update combine api client ([c1bb566](https://github.com/biosimulations/Biosimulations/commit/c1bb56636207d966bbd41f2dcf697008c43a27cf))
* **dispatch-service:** create handler to extract metadata ([16b469a](https://github.com/biosimulations/Biosimulations/commit/16b469a479061d643ffa8541a5db55e95f687404))
* **dispatch-service:** process and create metadata for simulation runs ([99df19d](https://github.com/biosimulations/Biosimulations/commit/99df19d50c786992aad6d35da1443fa3b8126c98))
* **exceptions:** change database errors to return 500 errors instead of 400 ([7390b7f](https://github.com/biosimulations/Biosimulations/commit/7390b7f3e35c6c1f0d1dfb8f145cf4ad8e5545c9))

## [3.13.0](https://github.com/biosimulations/Biosimulations/compare/v3.12.0...v3.13.0) (2021-08-24)


### Features

* **ontology:** updated to KiSAO 2.27 ([09ede72](https://github.com/biosimulations/Biosimulations/commit/09ede7208ce04542e5243ec96b4480471e98f8eb))

## [3.12.0](https://github.com/biosimulations/Biosimulations/compare/v3.11.0...v3.12.0) (2021-08-23)


### Bug Fixes

* **dispatch-api:** fields paramter is optional ([c3863e3](https://github.com/biosimulations/Biosimulations/commit/c3863e3c92709fb094ea6a9712aacbe59cdd412b))


### Features

* **dispatch-api:** implement metadata endpoints ([9d067e9](https://github.com/biosimulations/Biosimulations/commit/9d067e983cd625a8d706bc1cb3cfa2033bdabf62))

## [3.11.0](https://github.com/biosimulations/Biosimulations/compare/v3.10.0...v3.11.0) (2021-08-22)


### Features

* **combine-service:** enabled NEURON, NetPyNe for simulation ([19661df](https://github.com/biosimulations/Biosimulations/commit/19661df2b96697934d2f9ca3b9949cac1570554e))

## [3.10.0](https://github.com/biosimulations/Biosimulations/compare/v3.9.0...v3.10.0) (2021-08-21)


### Features

* **combine-service:** added endpoint for low latency simulation ([b44f5e3](https://github.com/biosimulations/Biosimulations/commit/b44f5e30df6b0e23b40ed302fa6a7a61a6969c09))
* **combine-service:** added endpoint for low-latency simulations ([931f3c7](https://github.com/biosimulations/Biosimulations/commit/931f3c7d25b329f55578cb0e7ba5e59ce7c91858))
* **combine-service:** added options to export simulation results in HDF5, zip formats ([c0cd699](https://github.com/biosimulations/Biosimulations/commit/c0cd6998dc850c5afbcee00fa2afec7390ffd235))
* **combine-service:** added simulator name to get simulaton tools endpoint ([1f708dd](https://github.com/biosimulations/Biosimulations/commit/1f708dd70e138a2c7d662b8689ef73566180e3ba))
* **combine-service:** added test to verify simulator APIs ([0bcaff3](https://github.com/biosimulations/Biosimulations/commit/0bcaff3f11fbc9e85285ef709432219d80fb1ef1))
* **combine-service:** pre-compiled Python code for faster initial calls ([f368de9](https://github.com/biosimulations/Biosimulations/commit/f368de9d973274c0e218e24075665a255d698e7f))
* **combine-service:** updated AMICI, GillesPy2, libSBMLSim ([fe8df2a](https://github.com/biosimulations/Biosimulations/commit/fe8df2a7b1da8c955210be78432565f60080d5bd))
* **combine-service:** updated to KiSAO 2.26, BioSimalators-utils 0.1.105 ([05b9dc5](https://github.com/biosimulations/Biosimulations/commit/05b9dc5ebe393c1e4c4143a0125867a3a79a58ea))
* **combine-service:** updated to KiSAO 2.26, BioSimulators-utils 0.1.105 ([74d0b36](https://github.com/biosimulations/Biosimulations/commit/74d0b369a933ec4b3f7b46c4028addfb794ef90d))
* **dispatch:** added simulator names to simulation tools menu in run form ([5170810](https://github.com/biosimulations/Biosimulations/commit/51708108ce772851690a4a3edff556513feb9368))
* **ontology:** updated to KiSAO 2.26 ([f25f243](https://github.com/biosimulations/Biosimulations/commit/f25f243bf85af086902badfd22b25547d974b9fb))
* **simulators:** added documentation for Python API conventions ([60669be](https://github.com/biosimulations/Biosimulations/commit/60669be5f166ac30a1af137240fbc049c13db331))

## [3.9.0](https://github.com/biosimulations/Biosimulations/compare/v3.8.0...v3.9.0) (2021-08-19)


### Features

* **dispatch:** added example run for MASSpy ([a589891](https://github.com/biosimulations/Biosimulations/commit/a589891a9367b58be7317f6b8f8da8545c2c44a7))
* **dispatch,ontology:** started to add MASS, RBA formats ([43a6153](https://github.com/biosimulations/Biosimulations/commit/43a615325f9695700ac2c9b68b2e124b0b03e3f9))
* **ontology:** updated to KiSAO 2.25 ([3fb5c54](https://github.com/biosimulations/Biosimulations/commit/3fb5c54af97d5631228c0e6456390377a867de5b))

## [3.8.0](https://github.com/biosimulations/Biosimulations/compare/v3.7.0...v3.8.0) (2021-08-18)


### Features

* **ontology:** updated to kisao 2.26 ([0f1f31a](https://github.com/biosimulations/Biosimulations/commit/0f1f31ae0383af38f9e7cd06aa28f022b7d6df07))

## [3.7.0](https://github.com/biosimulations/Biosimulations/compare/v3.6.0...v3.7.0) (2021-08-18)


### Bug Fixes

* update angular/cdk ([c7a18c1](https://github.com/biosimulations/Biosimulations/commit/c7a18c16f0c71fb0057a8cd331bbdd28abddb74b))
* **deps:** update dependency @openapi-contrib/openapi-schema-to-json-schema to v3.1.1 ([#2799](https://github.com/biosimulations/Biosimulations/issues/2799)) ([5d6453f](https://github.com/biosimulations/Biosimulations/commit/5d6453f1fc211b61ba7a0b40c9e18e7877bb874d))
* **deps:** update dependency @sendgrid/mail to v7.4.6 ([#2800](https://github.com/biosimulations/Biosimulations/issues/2800)) ([1e83398](https://github.com/biosimulations/Biosimulations/commit/1e83398a2b75f5af5b44506603f906f3c6e00b9d))
* **deps:** update dependency @typegoose/typegoose to v7.6.3 ([#2803](https://github.com/biosimulations/Biosimulations/issues/2803)) ([1a02d14](https://github.com/biosimulations/Biosimulations/commit/1a02d14610ec792e967f71cd4eec35c83af31d74))
* **deps:** update dependency auth0 to v2.36.1 ([#2823](https://github.com/biosimulations/Biosimulations/issues/2823)) ([fa11231](https://github.com/biosimulations/Biosimulations/commit/fa11231f161bb778e67f52287d707756b5315f20))
* **deps:** update dependency cache-manager to v3.4.4 ([#2804](https://github.com/biosimulations/Biosimulations/issues/2804)) ([3c62abb](https://github.com/biosimulations/Biosimulations/commit/3c62abb0576ede01ba5d9eec1ca0daf56a8b5a4e))
* **deps:** update nest ([#2807](https://github.com/biosimulations/Biosimulations/issues/2807)) ([875673f](https://github.com/biosimulations/Biosimulations/commit/875673f37fe2e4d0cc3f33497043f27899d03fd9))
* **dispatch-service:** added handling for case when no environment variables need to be set ([b51418d](https://github.com/biosimulations/Biosimulations/commit/b51418db31f38c534e774603070943d2abd0901e))
* **platform:** fix import of simulationrun metadata ([fb85650](https://github.com/biosimulations/Biosimulations/commit/fb85650554373a20e48be805d1068f0885bee7c0))


### Features

* **datamodel:** add common api query parameters ([c8bace5](https://github.com/biosimulations/Biosimulations/commit/c8bace5a78d6e991f0590a403927ab28757895f9))
* **dispatch:** added support for passing environment variables to simulators ([107221a](https://github.com/biosimulations/Biosimulations/commit/107221acfef0d611dba986f1da1e1782c1472d91))
* **dispatch-api:** add ability to get sparse simulationRuns ([5570ebb](https://github.com/biosimulations/Biosimulations/commit/5570ebb053df9d19c1fcb8f838564b1107d4c33f))
* **dispatch-api:** add endpoint for metadata ([24ffa40](https://github.com/biosimulations/Biosimulations/commit/24ffa40e96f42054ae4cee1f3eb62c6544ff08fe))
* **dispatch-api:** add tags, add ontology endpoint ([d056b4a](https://github.com/biosimulations/Biosimulations/commit/d056b4a74f0d19454e0eff0c515aea2eba5a0cd5))
* **dispatch-api:** create api model for metadata ([#2815](https://github.com/biosimulations/Biosimulations/issues/2815)) ([d55f7a5](https://github.com/biosimulations/Biosimulations/commit/d55f7a5d29e1f3cbb23ab472f84a2d8b961af843))
* **exceptions:** add better handling of validation errors ([e4e1986](https://github.com/biosimulations/Biosimulations/commit/e4e198623ae0d3f86460da18040faf9f28c9ee9a))

## [3.6.0](https://github.com/biosimulations/Biosimulations/compare/v3.5.0...v3.6.0) (2021-08-11)


### Bug Fixes

* **dispatch:** add logging when catching error ([85fab6c](https://github.com/biosimulations/Biosimulations/commit/85fab6cd89789cbf7357569194404e824bb60865))


### Features

* **dispatch:** added example runs for represillator model with SBML ([3400fa1](https://github.com/biosimulations/Biosimulations/commit/3400fa13eed1178ba74f76111b8ec83c995580f9))
* **dispatch:** added example simulation run for represillator model with OpenCOR ([ebffbae](https://github.com/biosimulations/Biosimulations/commit/ebffbae270595afd8d479d0a5ab6e90623b2323a))
* **dispatch:** added example simulation runs with visuaulizations using SBGN PD maps ([3c371e0](https://github.com/biosimulations/Biosimulations/commit/3c371e00758ab7a56b54afb5997da485fa3d071c))

## [3.5.0](https://github.com/biosimulations/Biosimulations/compare/v3.4.1...v3.5.0) (2021-08-09)


### Bug Fixes

* **deps:** pin dependency zone.js to 0.11.4 ([51bd9c7](https://github.com/biosimulations/Biosimulations/commit/51bd9c77245528ba8903d1644ddf5985899aa803))
* **dispatch:** added timeout for loading similar algorithms ([8c7a725](https://github.com/biosimulations/Biosimulations/commit/8c7a7256e13c9d6fd6c88fd4b5c89553d14906ef))
* **dispatch:** corrected display metadata loading indicator ([346aaeb](https://github.com/biosimulations/Biosimulations/commit/346aaeb809007d25b27e9a970f22eb7a3a715069))
* **dispatch:** fixed display of visualization loading indicator ([59fcd98](https://github.com/biosimulations/Biosimulations/commit/59fcd98686da9f8c6d22af7048b178d1b3362611))


### Features

* **dispatch:** added buttons to publish projects ([33c199a](https://github.com/biosimulations/Biosimulations/commit/33c199ae9e9ccdcc3ef41778f83c7caa7290d65e))
* **dispatch:** added display of errors with metadata of COMBINE archives ([41a0146](https://github.com/biosimulations/Biosimulations/commit/41a0146c8ac7421f13f0a9fe8aa61d2018f38ff0))
* **dispatch:** added support for XPP ([3d108ad](https://github.com/biosimulations/Biosimulations/commit/3d108adb5e468d408b4fbad061eab48f84861ea0))
* **dispatch:** improving capabilities when COMBINE service is down ([dc9ca09](https://github.com/biosimulations/Biosimulations/commit/dc9ca095ff96405e46a13d5f7ce56b95645308e7))
* **dispatch:** simplified designing 2D line/scatter plots with multiple curves ([beab68e](https://github.com/biosimulations/Biosimulations/commit/beab68e194e637940c57637520f292d0262fa328))

## [3.4.1](https://github.com/biosimulations/Biosimulations/compare/v3.4.0...v3.4.1) (2021-07-30)


### Bug Fixes

* **deps:** pin dependency @ngbmodule/material-carousel to 0.7.1 ([0d1acde](https://github.com/biosimulations/Biosimulations/commit/0d1acde594b0bc455c50209603684b8da7f66a02))
* **dispatch-api:** send correct message when simulation status changes ([a3c9c62](https://github.com/biosimulations/Biosimulations/commit/a3c9c6235102a33dafbc8414d0d5535c1a641f2f)), closes [#2739](https://github.com/biosimulations/Biosimulations/issues/2739)

## [3.4.0](https://github.com/biosimulations/Biosimulations/compare/v3.3.0...v3.4.0) (2021-07-29)


### Bug Fixes

* **dispatch:** corrected processing of metadata while status is pinging ([e29fb0e](https://github.com/biosimulations/Biosimulations/commit/e29fb0e54343da3c43db5f949dc48e83842845ca))


### Features

* **dispatch:** added example simulation run for activity flow diagram ([f534c33](https://github.com/biosimulations/Biosimulations/commit/f534c33835274eb113c9f62209d13810cb55f778))
* **dispatch:** expanded support for connecting SED-ML to Vega ([439bbeb](https://github.com/biosimulations/Biosimulations/commit/439bbebcfcc562304018711b9e0e485cba099eb1))
* **ontology:** updated SBO for additional framework terms ([5eaa097](https://github.com/biosimulations/Biosimulations/commit/5eaa097753353e9134a81b92554f3ca7efd8335e))

## [3.3.0](https://github.com/biosimulations/Biosimulations/compare/v3.2.0...v3.3.0) (2021-07-23)


### Bug Fixes

* **dispatch:** hiding figures/tables section when there are no figures/tables ([056caf8](https://github.com/biosimulations/Biosimulations/commit/056caf8d6f7ec6c7b0f6ad2116588e8f6dea751d))


### Features

* **dispatch,simulators:** added documentation about generating data visualizations ([0066522](https://github.com/biosimulations/Biosimulations/commit/00665225875e75a62a93dd93fd545f8e823f9ecc))
* making creation data metadata optional ([7812d65](https://github.com/biosimulations/Biosimulations/commit/7812d654f1ddfed9f9c2ea00b63ae31c3a537942))

## [3.2.0](https://github.com/biosimulations/Biosimulations/compare/v3.1.0...v3.2.0) (2021-07-22)


### Bug Fixes

* **dispatch-api:** change path from 'run' to 'runs' ([ead8d80](https://github.com/biosimulations/Biosimulations/commit/ead8d807ffe48a2cb50e54d88a65e880d00b6a70))


### Features

* **dispatch:** improved Vega error handling ([56a1e0c](https://github.com/biosimulations/Biosimulations/commit/56a1e0ce1f7147130bd18651dcac9d0b6953bb09))
* **dispatch:** updated example runs for new vis and metadata ([1683451](https://github.com/biosimulations/Biosimulations/commit/168345199fb240e098f168211609973076251a0b))
* **platform,platform-api:** platform gets projects from api ([f0b010d](https://github.com/biosimulations/Biosimulations/commit/f0b010d68b592765acb172c27a1b527ca4d9d157))
* **simulators:** added documentation about recommendation to use Identifiers.org URIs ([5445b08](https://github.com/biosimulations/Biosimulations/commit/5445b08a46678b4f73ea25fec53b38c9fdc6de4d))

## [3.1.0](https://github.com/biosimulations/Biosimulations/compare/v3.0.2...v3.1.0) (2021-07-19)


### Bug Fixes

* **deps:** pin dependencies ([5c762b6](https://github.com/biosimulations/Biosimulations/commit/5c762b64117e22863f8edf96c0d358256536b765))
* **deps:** update nest monorepo ([790aa52](https://github.com/biosimulations/Biosimulations/commit/790aa52b225b7d1eca6312da65fa0ff6b3d6fb9c))
* **dispatch:** fixed handling of query arguments to run simulation route ([546c49b](https://github.com/biosimulations/Biosimulations/commit/546c49bc344dfade3fd677998a665066e35ea2ce))
* **simulators:** reenabling display of simulator validation test results ([902bc93](https://github.com/biosimulations/Biosimulations/commit/902bc9374875572f4072ffde72d369a60f6f31c1)), closes [#2696](https://github.com/biosimulations/Biosimulations/issues/2696)


### Features

* **dispatch:** added 1D histogram plot ([379e0a7](https://github.com/biosimulations/Biosimulations/commit/379e0a7d74ad74b1cdc35e7babba9be7b1155c01))
* **dispatch:** added 2D heatmap data visualization ([147b6ad](https://github.com/biosimulations/Biosimulations/commit/147b6ad4ae52e77b4c02eac6cfc5b17dacf8a5c1))
* **dispatch:** added ability to add files to COMBINE archives ([c05631e](https://github.com/biosimulations/Biosimulations/commit/c05631e52fbcba6c72ccb24380b2d66167790218))
* **dispatch:** added exporting user-configured histogram viz to Vega ([6342e44](https://github.com/biosimulations/Biosimulations/commit/6342e448059fcc100d37b7246206cfb01b1bd8e8))
* **dispatch:** added Vega export for 2D heatmap, improved visualization form validation ([784917f](https://github.com/biosimulations/Biosimulations/commit/784917f3172e2416fe3e959218ec27e452fe4a79))
* **dispatch:** added Vega export for 2d line/scatter plot ([8f2cff0](https://github.com/biosimulations/Biosimulations/commit/8f2cff046bc3db7fdab6a8690a780da1d3ae7867))
* **dispatch:** improved plotting ([852545b](https://github.com/biosimulations/Biosimulations/commit/852545b75fb4dd1bca2b0594ca519f1cab9a111d))
* **dispatch:** linked Vega signals to attributes of SED-ML simulations ([88a68c3](https://github.com/biosimulations/Biosimulations/commit/88a68c332c9f630bc9c5cffa02b9c9e2e3f0058a))
* **dispatch:** linking published figures/tables to displayed visualizations ([44c810a](https://github.com/biosimulations/Biosimulations/commit/44c810acdc8ed31173387bd5521dbc03c093008a))
* **platform:** implement viewing a project ([0ad9af3](https://github.com/biosimulations/Biosimulations/commit/0ad9af33635272cfc29e4cd9b3d2d71cdb03dbe4))
* **platform-api:** add skeleton implementation ([0758052](https://github.com/biosimulations/Biosimulations/commit/0758052e88b185d72e3cedc46256377a8e3d9753))

## [3.0.2](https://github.com/biosimulations/Biosimulations/compare/v3.0.1...v3.0.2) (2021-07-13)


### Bug Fixes

* **mail-service,dispatch-service:** fix nats-server connection ([9143d0c](https://github.com/biosimulations/Biosimulations/commit/9143d0cabba6a22d44398cd4eec1cb1a5033e2f9))

## [3.0.1](https://github.com/biosimulations/Biosimulations/compare/v3.0.0...v3.0.1) (2021-07-13)


### Bug Fixes

* try new server options ([1851742](https://github.com/biosimulations/Biosimulations/commit/1851742fdb222a6330a6b2f52b814ee7d3273c5a))

## [3.0.0](https://github.com/biosimulations/Biosimulations/compare/v2.5.2...v3.0.0) (2021-07-13)


### Bug Fixes

* **mail-service,dispatch-service:** fix import of http module ([805e48f](https://github.com/biosimulations/Biosimulations/commit/805e48f34dec1b50308c20affa786fac1ae646f5))


### Features

* **simulators-api:** add a query argument to include the results of the validation ([710be08](https://github.com/biosimulations/Biosimulations/commit/710be085ec732d851aa89e78773c4ba12e7e682e)), closes [#2668](https://github.com/biosimulations/Biosimulations/issues/2668)


### BREAKING CHANGES

* **simulators-api:** validation data is no longer returned by default. A Query argument is needed to include the validation information

## [2.5.2](https://github.com/biosimulations/Biosimulations/compare/v2.5.1...v2.5.2) (2021-07-13)


### Bug Fixes

* type and build fixes ([6812bd0](https://github.com/biosimulations/Biosimulations/commit/6812bd0de12c3716ac8928154ed73aa92953dc40))
* **dispatch-api:** fix error with parsing outputIds ([9fac99f](https://github.com/biosimulations/Biosimulations/commit/9fac99f68094919de23d632795bf131b9fb8a1ef)), closes [#2683](https://github.com/biosimulations/Biosimulations/issues/2683)

## [2.5.1](https://github.com/biosimulations/Biosimulations/compare/v2.5.0...v2.5.1) (2021-07-09)


### Bug Fixes

* **simulators-api:** Allow for date based versions ([0c8fb8d](https://github.com/biosimulations/Biosimulations/commit/0c8fb8d00a36675f652a665d9c279709e798c212)), closes [#2681](https://github.com/biosimulations/Biosimulations/issues/2681)

## [2.5.0](https://github.com/biosimulations/Biosimulations/compare/v2.4.0...v2.5.0) (2021-07-09)


### Bug Fixes

* **dispatch:** downloading created COMBINE/OMEX archives ([03895bf](https://github.com/biosimulations/Biosimulations/commit/03895bf82abfa0b6d9c7c7db186a31d29e726b49))
* correcting size of form fields ([ae69630](https://github.com/biosimulations/Biosimulations/commit/ae69630cf3e45abfbcaaec17787956e7189e5e53))
* **shared-exceptions:** include error metadata in the "meta" output ([2be0178](https://github.com/biosimulations/Biosimulations/commit/2be0178af9003156ad25fa22e8c7fe51457c9556))


### Features

* **combine-service:** adding support for creating steady-state analyses of logical models ([a9e6667](https://github.com/biosimulations/Biosimulations/commit/a9e6667034c2b35d4379dff72a2d7cefe4d4f4d8))
* **combine-service:** updating to biosimulators-utils 0.1.93 ([ca0a21e](https://github.com/biosimulations/Biosimulations/commit/ca0a21e33d7c2a54f8bd6d9aa9d8c6943da955b2))
* **shared-exceptions:** add validation pipe error factory ([35edb4d](https://github.com/biosimulations/Biosimulations/commit/35edb4d4f73d82e4bbb17bbd701d13fc580093af))

## [2.4.0](https://github.com/biosimulations/Biosimulations/compare/v2.3.0...v2.4.0) (2021-07-08)


### Bug Fixes

* **combine-service:** add protocol to server in API spec ([bcf4119](https://github.com/biosimulations/Biosimulations/commit/bcf41192f15894f7993239f16b14f415c1a85910))
* **combine-service:** debugged specifications ([f8b9420](https://github.com/biosimulations/Biosimulations/commit/f8b9420fadbcf10c367f26c6b77bf736dc5463f3))
* **combine-service:** fix api spec ([ebab457](https://github.com/biosimulations/Biosimulations/commit/ebab457a262b594744e81cfdfe4bfdd9af45a4c2))
* **platform:** enable strict template checking, fix type errors ([9facdb1](https://github.com/biosimulations/Biosimulations/commit/9facdb19c1cacae3d22d4bd952c0f1f0cbaf8035)), closes [#2185](https://github.com/biosimulations/Biosimulations/issues/2185)


### Features

* **combine-service:** add combine-service api client library ([bfa25b8](https://github.com/biosimulations/Biosimulations/commit/bfa25b8f0319de62d3c6a5902597d51c68b8eb96))
* **dispatch:** added example simulation run for GINsim ([3e639b3](https://github.com/biosimulations/Biosimulations/commit/3e639b30c49e3149d966518d5a82908a3905e831))
* **dispatch:** added example simulation runs for GINsim, LibSBMLSim ([06ad0a4](https://github.com/biosimulations/Biosimulations/commit/06ad0a4be78dc595305492e975ba856722d27b27))
* **dispatch:** adding support for GINML, ZGINML to COMBINE archive creation and execution ([9f949e5](https://github.com/biosimulations/Biosimulations/commit/9f949e561256288f2851468a83c96160cc14f7fe))
* **dispatch,ontology:** add terms for GINsim format ([22d8a7b](https://github.com/biosimulations/Biosimulations/commit/22d8a7b2daec93c086bacc3c61a279dc85481cfd))
* **ontology:** updating to KiSAO 2.19 with terms for logical modeling ([c47e63b](https://github.com/biosimulations/Biosimulations/commit/c47e63b7767eecde90c033d8c11ef55b89678d4a))
* **ontology,combine-service:** update to KiSAO 2.20 ([fadf3da](https://github.com/biosimulations/Biosimulations/commit/fadf3da7c0e714267a89ac903acf516be5f00533))

## [2.3.0](https://github.com/biosimulations/Biosimulations/compare/v2.2.1...v2.3.0) (2021-07-02)

### Features

- **dispatch:** Added tab to simulation run page to display metadata about the simulation project ([#2667](https://github.com/biosimulations/Biosimulations/issues/2667)) ([dde87fa](https://github.com/biosimulations/Biosimulations/commit/dde87faae5e558c3bbe86f6f17467ae747da55d8)), closes [#2661](https://github.com/biosimulations/Biosimulations/issues/2661)

## [2.2.1](https://github.com/biosimulations/Biosimulations/compare/v2.2.0...v2.2.1) (2021-07-01)

### Bug Fixes

- **dispatch:** fix example simulation runs ([60d91c1](https://github.com/biosimulations/Biosimulations/commit/60d91c1bb70e6ae08274a9380143baa19fa51043)), closes [#2653](https://github.com/biosimulations/Biosimulations/issues/2653)
- **simulators-api:** fix getting latest version ([4594c96](https://github.com/biosimulations/Biosimulations/commit/4594c96b53859e03960458cd001cf8614d64f64c)), closes [#2664](https://github.com/biosimulations/Biosimulations/issues/2664)

## [2.2.0](https://github.com/biosimulations/Biosimulations/compare/v2.1.0...v2.2.0) (2021-06-30)

### Bug Fixes

- **dispatch:** correct integration between simulation results and SED plots ([0bab60f](https://github.com/biosimulations/Biosimulations/commit/0bab60fe06cc52d55a670d8957e385dc7f247854))
- **dispatch:** download file instead of redirect ([cd2840d](https://github.com/biosimulations/Biosimulations/commit/cd2840d98d84f13eab34cea479a09da23187fe14)), closes [#2435](https://github.com/biosimulations/Biosimulations/issues/2435)
- **dispatch:** use correct api to get simulator info ([1e66f1f](https://github.com/biosimulations/Biosimulations/commit/1e66f1f85f4436987ca034c3cdafad9536c12b9e))

### Features

- add some shared endpoints ([567e4c2](https://github.com/biosimulations/Biosimulations/commit/567e4c27de05655d3b78b441e84231977afd234b))
- **auth-client:** Cache tokens locally ([f53c9f8](https://github.com/biosimulations/Biosimulations/commit/f53c9f8d4c9c3e2bed497ec85c4c53d774af9fb1)), closes [#2503](https://github.com/biosimulations/Biosimulations/issues/2503)
- **auth-common:** add util functions ([e0ac842](https://github.com/biosimulations/Biosimulations/commit/e0ac842518af8e6909493cb1b2b774a56faf6b17))

### Performance Improvements

- **dispatch-service:** use /local as the singularity cache/working directory ([c63b58c](https://github.com/biosimulations/Biosimulations/commit/c63b58c35a0c3da71910523a3baf0f445f5e493a))

## [2.1.0](https://github.com/biosimulations/Biosimulations/compare/v2.0.0...v2.1.0) (2021-06-18)

### Bug Fixes

- **dispatch-service:** remove check for process flag ([f7f88cc](https://github.com/biosimulations/Biosimulations/commit/f7f88cce2fbc54df13e34ef5212f1491036ec8b5)), closes [#2577](https://github.com/biosimulations/Biosimulations/issues/2577)

### Features

- **dispatch-api, dispatch-service:** add status reason to datamodel ([ca9bcb6](https://github.com/biosimulations/Biosimulations/commit/ca9bcb6c7d7ffcb0328ef679d5a82801995add45)), closes [#2441](https://github.com/biosimulations/Biosimulations/issues/2441)

## [2.0.0](https://github.com/biosimulations/Biosimulations/compare/v1.0.0...v2.0.0) (2021-06-17)

### Bug Fixes

- **dispatch-api:** bind class to this variable in map ([b4bb3ca](https://github.com/biosimulations/Biosimulations/commit/b4bb3ca27cd52d27abe68dcaa524a158a1a73507))
- dispatch frontend uses the updated api parameter ([#2636](https://github.com/biosimulations/Biosimulations/issues/2636)) ([a13779c](https://github.com/biosimulations/Biosimulations/commit/a13779cdc320d58c595f85399ca4d7747d603657)), closes [#2635](https://github.com/biosimulations/Biosimulations/issues/2635)

### Features

- **dispatch-api, dispatch-service:** Use HSDS to get simulation run data ([33b8030](https://github.com/biosimulations/Biosimulations/commit/33b8030e60fcbd2eb693e2a962620cf42855b4e4)), closes [#2533](https://github.com/biosimulations/Biosimulations/issues/2533) [#2442](https://github.com/biosimulations/Biosimulations/issues/2442) [#2440](https://github.com/biosimulations/Biosimulations/issues/2440) [#2369](https://github.com/biosimulations/Biosimulations/issues/2369) [#2069](https://github.com/biosimulations/Biosimulations/issues/2069)

### BREAKING CHANGES

- **dispatch-api, dispatch-service:** Dispatch API no longer has endpoints for creating or updating "Result" objects.
  The output of the results endpoints are updated to include information about type and shape of the data.
  The parameter "sparse" has been changed to "includeData".
  The datamodel for results has been adjusted to include all outputs, not just reports. "reports" has been renamed to "outputs"

# 1.0.0 (2021-06-16)

This is an arbitrary starting point for tracking changes and versioning. It should not be considered as the "first release".

### Bug Fixes

- bash script ([866b58a](https://github.com/biosimulations/Biosimulations/commit/866b58a244d3483fa6afa2ae0e8383e234920ba6))
- add check for large files downloading ([ca10aa5](https://github.com/biosimulations/Biosimulations/commit/ca10aa5f2c44fafcff8471c0da0bbd9db2a655ed)), closes [#2536](https://github.com/biosimulations/Biosimulations/issues/2536)
- bring inline with datamodel ([af53a54](https://github.com/biosimulations/Biosimulations/commit/af53a54a5e5bc91835e114e67e10fc69883c7f9b))
- change url to download results ([3d264d4](https://github.com/biosimulations/Biosimulations/commit/3d264d4abfdb85aa4cae99e1477cbf4666b8ba36)), closes [#2561](https://github.com/biosimulations/Biosimulations/issues/2561)
- check job status after completion ([c16649c](https://github.com/biosimulations/Biosimulations/commit/c16649c507c35f1e086d41fb496a573549f925ba))
- cleanup logs ([fbd330f](https://github.com/biosimulations/Biosimulations/commit/fbd330f2762f6192b7e269ffcd2e8ede93c5ad14))
- correct value for constant ([e2d3a68](https://github.com/biosimulations/Biosimulations/commit/e2d3a68a39f5e6ad3216daf109087a2b1c43f26b))
- fix error in reading port ([e1f6fb9](https://github.com/biosimulations/Biosimulations/commit/e1f6fb923a42283d1b42765b4d0376a146f406ef))
- fix logs and context buttons ([777e8e8](https://github.com/biosimulations/Biosimulations/commit/777e8e8f79f829b3762c7aa189a9d6184f4b24a1)), closes [#2543](https://github.com/biosimulations/Biosimulations/issues/2543) [#2540](https://github.com/biosimulations/Biosimulations/issues/2540)
- fix redis queue and port ([5f33a19](https://github.com/biosimulations/Biosimulations/commit/5f33a192203323e30d6badd4b6500cc056b3ef34))
- fix s3 key for downloading outputs ([f585a9a](https://github.com/biosimulations/Biosimulations/commit/f585a9a6295cbbf60e73c05d5ae908713d1ef5ee)), closes [#2622](https://github.com/biosimulations/Biosimulations/issues/2622)
- fix spelling of library ([a471e95](https://github.com/biosimulations/Biosimulations/commit/a471e95e6684ee093d036f343efbba50df327563))
- fix test ([6f236df](https://github.com/biosimulations/Biosimulations/commit/6f236df6b5186e44ccf9c459faa83698cf22d7ae))
- fix test ([6af0ca8](https://github.com/biosimulations/Biosimulations/commit/6af0ca8a0b9f2d557a0fd416475151261a46fb88))
- lint fix ([a26c24b](https://github.com/biosimulations/Biosimulations/commit/a26c24b17e9d1f72bdb53860b7f27b300030ec68))
- order of operations for creating results ([eac31e0](https://github.com/biosimulations/Biosimulations/commit/eac31e01bde327b1d3ff89f6a7cf7480e5d0c96d))
- propely set name and filetype of outputs ([951c239](https://github.com/biosimulations/Biosimulations/commit/951c239983ea93009565287a8ac9bfd3deae8052))
- Remove bad library import ([ecc86fa](https://github.com/biosimulations/Biosimulations/commit/ecc86fa6d9abf59b466ea02d25d62c1119d07de8)), closes [#2420](https://github.com/biosimulations/Biosimulations/issues/2420)
- remove xdg runtime directory ([f5ec15b](https://github.com/biosimulations/Biosimulations/commit/f5ec15bd726ab4afa01b0c2be4688217d4d89198))
- resolve build errors ([6691ebe](https://github.com/biosimulations/Biosimulations/commit/6691ebedbda107862cbf731cb891044c426e5fc9))
- typo in return statement ([1f6c4fc](https://github.com/biosimulations/Biosimulations/commit/1f6c4fc0bda0780170530789e27e4fab4233d2e3))
- update default stoage URL ([f9b0d75](https://github.com/biosimulations/Biosimulations/commit/f9b0d75b3c4370653d1c0596157cc381fa6573f0))
- update logs ([818a0c3](https://github.com/biosimulations/Biosimulations/commit/818a0c347529c42d697ac972c15e17c09b5e0372))
- update sbatch memoy amount ([b9026f9](https://github.com/biosimulations/Biosimulations/commit/b9026f96ff2e5b4876559b1b116c1d5cdebbfb8d))
- update sbatch script to use custom module ([0ef1c52](https://github.com/biosimulations/Biosimulations/commit/0ef1c52de4d6703032decaff9b2c8941175c70fb))
- use job status to determine completion ([adb12a0](https://github.com/biosimulations/Biosimulations/commit/adb12a0efbe07e82346ddada18ad93342d1cede5))
- **apps/frontend:** relative import ([3854f27](https://github.com/biosimulations/Biosimulations/commit/3854f272fdd21847e522cb03f25353b06a3c3028))
- **auth:** check for logged in before intercepting ([7c22a19](https://github.com/biosimulations/Biosimulations/commit/7c22a19a2a33cff63067d30dd19ff8bfe091a189))
- **auth:** Check for username correctly ([05996b3](https://github.com/biosimulations/Biosimulations/commit/05996b376ed7b08a4de974cf39029cb9956c1070))
- **disatch:** Fix open api schema ([ff15503](https://github.com/biosimulations/Biosimulations/commit/ff15503862e62412418f3144353e76a8cd877f7b))
- **dispatch:** Fix observable piping ([34e7086](https://github.com/biosimulations/Biosimulations/commit/34e7086261b20ffdf42b6edcf3112f972432ea79))
- **dispatch:** parsing results accounts for quotes ([d055a7b](https://github.com/biosimulations/Biosimulations/commit/d055a7ba328e1e4d4d5f8094a661377a8e5294f9)), closes [#2459](https://github.com/biosimulations/Biosimulations/issues/2459)
- **dispatch:** patch error handling ([d2d98e5](https://github.com/biosimulations/Biosimulations/commit/d2d98e57bd1632d8f289e2d9fd017443a653c8db))
- **dispatch:** remove bad environment variables ([3c31b7d](https://github.com/biosimulations/Biosimulations/commit/3c31b7de39d05b5bfaf640a1553a204267247eab)), closes [#2476](https://github.com/biosimulations/Biosimulations/issues/2476)
- **dispatch:** Simulation results not saved for some simulations and overall status doesn't reflect such errors ([#2428](https://github.com/biosimulations/Biosimulations/issues/2428)) ([acd2dff](https://github.com/biosimulations/Biosimulations/commit/acd2dff837834e6732f4b5074c433f90a9523d06)), closes [#2416](https://github.com/biosimulations/Biosimulations/issues/2416)
- use https for auth0 image ([19a4dcc](https://github.com/biosimulations/Biosimulations/commit/19a4dcc53d6572be6b60a9bb8a4d9db4bd89afc6))
- **forms:** connect taxon form properly ([a4088f7](https://github.com/biosimulations/Biosimulations/commit/a4088f79bf3a6c06cfc0dfbb48936d820428c378))
- **forms:** fix reference form implementation ([5b3eba4](https://github.com/biosimulations/Biosimulations/commit/5b3eba498a464fe77f7bbad3ae3365db8bb3e6bc))
- **forms:** fix some taxon form details ([8802ced](https://github.com/biosimulations/Biosimulations/commit/8802cedb84fc7c1afdc2b5cea244bc11c5a67399))
- **forms:** fix tags form ([b60b99c](https://github.com/biosimulations/Biosimulations/commit/b60b99cc6c8526a7403c585f8f79315c757972a5))
- **forms:** fix validation error with taxon form ([6cb2d43](https://github.com/biosimulations/Biosimulations/commit/6cb2d43e7edda5edff8310cc5bcffc5a95cedec8))
- **forms:** Form does not scroll over topbar ([70f3275](https://github.com/biosimulations/Biosimulations/commit/70f327598c2cacd278625dd78c16e93761974d54))
- **forms:** set file form disable properly ([c523d0e](https://github.com/biosimulations/Biosimulations/commit/c523d0e8230c14c3c64a0266107b7777851ae114))
- **gaurds:** Gaurd loads underconstruction pages by default ([9f7a810](https://github.com/biosimulations/Biosimulations/commit/9f7a810bd6aa000039fb7dd3e3e5b421f63118e3))
- **grid:** grid work with async resources ([d621cd3](https://github.com/biosimulations/Biosimulations/commit/d621cd3d11402aca909a005e376d214160f117b2))
- **interceptor:** API token, error handling ([78be137](https://github.com/biosimulations/Biosimulations/commit/78be1377afbf3877885744074fe8dccd0ffe1ca6))
- **interceptor:** Fixed a bug in the error handling of the interceptor ([1397c96](https://github.com/biosimulations/Biosimulations/commit/1397c966e4a83fd1e0b660daeef092617cecc106))
- **models:** Edit component calls subscribe ([50f15ac](https://github.com/biosimulations/Biosimulations/commit/50f15accb120fcdb136262541e386bacee630cb3))
- **navigation:** get username via async ([7415eba](https://github.com/biosimulations/Biosimulations/commit/7415ebae27841b6147931820869b8ef06a55cc89))
- **navigation:** have navigation work with async ([0da33ea](https://github.com/biosimulations/Biosimulations/commit/0da33eaf426d907115ba07636674261611e73a7c))
- **polyfill:** add back polyfill ([30e9c1d](https://github.com/biosimulations/Biosimulations/commit/30e9c1de6204f3bba371be0ad53c89b8a8939f1f))
- **resource service:** add query params to read ([d1d7d38](https://github.com/biosimulations/Biosimulations/commit/d1d7d38610582ed45769b5f3d404f9222cbd08cc))
- **resources:** Move more functionality into abstract class ([5e94323](https://github.com/biosimulations/Biosimulations/commit/5e94323cdf7405de6d80b69bb18a7217d02e9922))
- **serializer:** fix private public being flipped ([035241f](https://github.com/biosimulations/Biosimulations/commit/035241f949d25c38bf9767a85ee805df958d0395))
- **serializers:** improve serializers ([0d3f004](https://github.com/biosimulations/Biosimulations/commit/0d3f004175ece11a48c668959d9b5b3e217994ae))
- **serializers:** user serializer returns none for '' ([221ce41](https://github.com/biosimulations/Biosimulations/commit/221ce413570b04671da1e77581d441a423eb4647))
- **services:** Dont cast http reponse to resource ([b20e69d](https://github.com/biosimulations/Biosimulations/commit/b20e69d45dcdcaf69dedc41cc157dc9cfee14005))
- **simulations:** fix async view ([d2cad0f](https://github.com/biosimulations/Biosimulations/commit/d2cad0fd3230714bae9be8858389456bcb6ac9ca))
- **tests:** fix common test issues ([61c7621](https://github.com/biosimulations/Biosimulations/commit/61c76219e91b7b43675014ea7293dd20e3dc1cf2))
- **tests:** fix common test issues ([0c4e5f3](https://github.com/biosimulations/Biosimulations/commit/0c4e5f3f9fa58cc600a45b401a5a4ab3ee23c315))
- **user:** profile edit component properly creates user ([0bed682](https://github.com/biosimulations/Biosimulations/commit/0bed68288d6fcb1c885138ec83e0c7309f3415d3))
- **visualizations:** fix licence view ([f5d9973](https://github.com/biosimulations/Biosimulations/commit/f5d997394786baf2fac0fe4a0b6ea11c402b9881))

### Features

- add client library ([820647b](https://github.com/biosimulations/Biosimulations/commit/820647b4c13ddae30174232c1da0cd5f88990dc6))
- add config for queue ([c7ec4a1](https://github.com/biosimulations/Biosimulations/commit/c7ec4a1ed6b840e11968597f5c530bf6d0a15566))
- add hsds client module ([1730521](https://github.com/biosimulations/Biosimulations/commit/17305212950d8f27e3d577660de00828b90a5f2d))
- use org for getting latest simulator ([2f7c503](https://github.com/biosimulations/Biosimulations/commit/2f7c503b30f921f8e21919934774c13dec4113d8))
- **api:** add config to api ([931dcf5](https://github.com/biosimulations/Biosimulations/commit/931dcf518c7d57adaeb767bdc18e7fcee151be1a))
- **api:** add crud skeleton for routes ([64fce18](https://github.com/biosimulations/Biosimulations/commit/64fce18062fa43ea5c3063fe29a569d8b97f1d09))
- **api:** add open api spec generation ([659f8b4](https://github.com/biosimulations/Biosimulations/commit/659f8b4a51cfb201aea00fc7d5b2bfac6c7a9c14))
- **author form:** build out author form ([afe666a](https://github.com/biosimulations/Biosimulations/commit/afe666ae014b6463d53349913488d894f64f5aa6))
- **datamodel:** add gaurds ([2ffdf04](https://github.com/biosimulations/Biosimulations/commit/2ffdf04605f2ecb9dd41d7c7c681e5843ece23c2))
- **datamodel:** add more of the core datamodel ([9b28f83](https://github.com/biosimulations/Biosimulations/commit/9b28f834cc0a671e1628617d97acc6bae9068412))
- **datamodel:** add properties to format ([67b024c](https://github.com/biosimulations/Biosimulations/commit/67b024ccb1620937a711da7486bc328f9f698328)), closes [#462](https://github.com/biosimulations/Biosimulations/issues/462)
- **datamodel:** add url to ontology ([3c8a169](https://github.com/biosimulations/Biosimulations/commit/3c8a1699b279acd687b1fcd80a420bb509cd09fe))
- **datamodel:** redefine core objects as set of attriutes ([95ea0a7](https://github.com/biosimulations/Biosimulations/commit/95ea0a78c2c23617b8ed103de55e41fdf2059122))
- **datamodel:** redfine core resources as primary and secondary ([dd2fec4](https://github.com/biosimulations/Biosimulations/commit/dd2fec4ac6d0f8fce9bfa3d1232e95b1dd528986))
- **dispatch:** add a dispatch service ([152b3b0](https://github.com/biosimulations/Biosimulations/commit/152b3b06536428850af601187eaf7f243f45b4d6))
- **errors:** add default errors component ([a10b927](https://github.com/biosimulations/Biosimulations/commit/a10b92773d263c289e59bab261381fb49bf3f953))
- **errors:** Create 404 component ([e1c7d33](https://github.com/biosimulations/Biosimulations/commit/e1c7d334ccde6684284ba1cab0ce38b2f109c3e9))
- **errors:** Create errors module ([7e3cac4](https://github.com/biosimulations/Biosimulations/commit/7e3cac4e27b4104595f94ef269e478916336b168))
- **errors:** Create under construction component ([5d2fbe1](https://github.com/biosimulations/Biosimulations/commit/5d2fbe1405117f8acc52b7ae5cbc9c5b7d4f5974))
- **errors:** slight changes to underConstruction ([bb06012](https://github.com/biosimulations/Biosimulations/commit/bb06012172c2a86c014ad1176832fe68a3439428))
- **forms:** add a component for file inputs ([e843f87](https://github.com/biosimulations/Biosimulations/commit/e843f873114fde51c917d01f48f6ce3f3fe306ff))
- **forms:** add abstract array subform ([1d70e47](https://github.com/biosimulations/Biosimulations/commit/1d70e479618af355416e287f973b5b8bbdc3db67))
- **forms:** add authors and identifiers components ([d2e35e8](https://github.com/biosimulations/Biosimulations/commit/d2e35e8bca295aa5a536eea3fbabc02c9119c13a))
- **forms:** Add identifier form control ([dd80d37](https://github.com/biosimulations/Biosimulations/commit/dd80d37d640c5eae612d216f26e3472ab5607cd8))
- **forms:** add name form control ([8557559](https://github.com/biosimulations/Biosimulations/commit/8557559e2e2d85e0409e31aa9cb2faaff39165a3))
- **forms:** add required validator to fields ([ceaacb5](https://github.com/biosimulations/Biosimulations/commit/ceaacb53174505e27ab2e98da1d80db4d346a6bf))
- **forms:** add resource form skeleton ([e7a56d5](https://github.com/biosimulations/Biosimulations/commit/e7a56d5892c80205fba48809e1f89dc347f3abff))
- **forms:** add taxon form ([48a2e86](https://github.com/biosimulations/Biosimulations/commit/48a2e8605a3de64f75368fd8a4bd0bcd8b6d6785))
- **forms:** add to resource form implementation ([dc0b9ed](https://github.com/biosimulations/Biosimulations/commit/dc0b9ed74e1ec8435e6422d8623c0b5aa75344b4))
- **forms:** add username form ([7ba4e81](https://github.com/biosimulations/Biosimulations/commit/7ba4e814ca1b5155bd5ae067c9385dbe3e84a30c))
- **forms:** create a model form ([353d2ab](https://github.com/biosimulations/Biosimulations/commit/353d2ab55ca5c4d6d7a3539a55e09604888bc5ef))
- **forms:** Create descriptions form control ([87b0fdd](https://github.com/biosimulations/Biosimulations/commit/87b0fdd95967957dbcc2033e858df9f975d838cb))
- **forms:** Create edit-preview component ([5177f05](https://github.com/biosimulations/Biosimulations/commit/5177f05ee5d86c6f19ee735b9c181f939971150d))
- **forms:** create model format form ([af0ad55](https://github.com/biosimulations/Biosimulations/commit/af0ad55b1972146aa6313e89dd8359a7d5d6d604))
- **forms:** enable access form component ([741afa3](https://github.com/biosimulations/Biosimulations/commit/741afa323769ffc6209e1b7cd6056311e3625a65))
- **forms:** Enable Drag/Drop ([799147b](https://github.com/biosimulations/Biosimulations/commit/799147b2277a7cc06191a50eb3cd6f66cc589f13))
- **forms:** Finalize author form ([3f4f86c](https://github.com/biosimulations/Biosimulations/commit/3f4f86cfb56c7d3785362f384a1a3378310a9cfd))
- **forms:** Generalize single field controls ([0da8122](https://github.com/biosimulations/Biosimulations/commit/0da8122c33d69a36a0a8db810e60a40a9340404d))
- **forms:** implement licence form ([2d2d9b1](https://github.com/biosimulations/Biosimulations/commit/2d2d9b1340e30867321ce8fbcef86c1e7e3be8e4))
- **forms:** implement refrences form ([122fa33](https://github.com/biosimulations/Biosimulations/commit/122fa33703cbafcaef32a4665ad130ee8599cd29))
- **forms:** implement resource form ([85e89ea](https://github.com/biosimulations/Biosimulations/commit/85e89ea97e8990e23b9b9477e9017b3f289854a1))
- **forms:** implement tags form ([0b1e480](https://github.com/biosimulations/Biosimulations/commit/0b1e48018e8cc1f22620e775cbd836fb63e315bb))
- **forms:** improve disable handling ([88cace8](https://github.com/biosimulations/Biosimulations/commit/88cace866e81e97def3dbbbedeacaf4daf24f435))
- **forms:** styling ([32044ef](https://github.com/biosimulations/Biosimulations/commit/32044ef922ae05010aebe1c9e83016007e7b069d))
- **home:** Add sponsors section ([8f84a67](https://github.com/biosimulations/Biosimulations/commit/8f84a677f550655b17d425dc6c961276a9fcaca5))
- **logging:** Added logging ([18476ca](https://github.com/biosimulations/Biosimulations/commit/18476ca216a6f0ddd5274d51073e7733c2be0c15))
- **login:** add styling ([642a71c](https://github.com/biosimulations/Biosimulations/commit/642a71cb269498b8272dd0a0b0e85c2f395fa170))
- **login:** redirect works ([eb05030](https://github.com/biosimulations/Biosimulations/commit/eb05030d658ff4570eac67c4fb08a91e38d70950))
- **mateiral:** add a material topbar ([4325dcc](https://github.com/biosimulations/Biosimulations/commit/4325dcc4fcb05cfd46a900b2df88bd33ce0ec32f))
- **models:** add query-options model definition ([d628289](https://github.com/biosimulations/Biosimulations/commit/d6282893f47393846af613e1cc4173c0f3b550ef))
- **projects:** add view project ([dab3c59](https://github.com/biosimulations/Biosimulations/commit/dab3c5917c4905c617e1fed78a50854aabc80d37))
- **pwa:** Add pwa capabilities ([27b9050](https://github.com/biosimulations/Biosimulations/commit/27b90508a31e840acc16fbccae809b2b288fde3b))
- **resources:** Resources now have owner embedded ([f03c30b](https://github.com/biosimulations/Biosimulations/commit/f03c30bf7ef7ad5096b27aaf431730c0e402c6ac))
- **serializers:** serializers read files properly ([c0eb5c1](https://github.com/biosimulations/Biosimulations/commit/c0eb5c15dae85aa9e164643e04c957359afca6fa))
- **service:** Breadcrumb service generates breadcrums ([68efc1a](https://github.com/biosimulations/Biosimulations/commit/68efc1a43e0ff52250b16edc73b4351c0b6f647c))
- **services:** add config service ([98ed473](https://github.com/biosimulations/Biosimulations/commit/98ed473bfa04ba40358ef53b376932a660abc286))
- **services:** add file service ([0538686](https://github.com/biosimulations/Biosimulations/commit/053868605725bc0e5fad0944c53502737e46b950))
- **services:** add test organism to metadataservice ([91f2503](https://github.com/biosimulations/Biosimulations/commit/91f25039a9b8d39781bd7b7a2137362de18c7358))
- **services:** Resource services return new model ([0c32ba5](https://github.com/biosimulations/Biosimulations/commit/0c32ba582802824e1dd4e9f64f332e145ead1bee))
- **shared:** add an authentication library ([904cd59](https://github.com/biosimulations/Biosimulations/commit/904cd59b320a33ec8447980e79c12bbcb3121a29))
- **shared:** add fields to remote file ([b85ff9d](https://github.com/biosimulations/Biosimulations/commit/b85ff9d7061a1e62d191b017078e680d5429c7b9))
- **shared:** add more model serializing ([b06129d](https://github.com/biosimulations/Biosimulations/commit/b06129d1b0b7e2b4499c94fb34490b38a15345a3))
- **visualization:** More flexible 2D visualization to better match needs of SED-ML L1V3 and BioModels ([31c96de](https://github.com/biosimulations/Biosimulations/commit/31c96de50f6c1323b574e00f400e5e13ade93060))
- add construction gaurd in prodction mode ([84708a2](https://github.com/biosimulations/Biosimulations/commit/84708a21f656318260e6162518fef84997051916))
- create a debugger component ([fd47c37](https://github.com/biosimulations/Biosimulations/commit/fd47c374f70d4dc28646cec31f995075c778d7a3))
- **shared:** add under constrcution gaurd ([d931f6d](https://github.com/biosimulations/Biosimulations/commit/d931f6de41fc1adfa558355664b651ce4a160c6d))
- **shared:** Remote file has a method to create from File ([f3360cf](https://github.com/biosimulations/Biosimulations/commit/f3360cfb6cfbb733d46501327986e7c9aec565d2))
- **simulations:** Simulations now have embedded models ([4dd512b](https://github.com/biosimulations/Biosimulations/commit/4dd512bc1944a6a2b2e8b00a8d5f48bd61939e9f))
- **users:** Add different snackbars ([056ffe2](https://github.com/biosimulations/Biosimulations/commit/056ffe2e4f26b4f43835a9cb7978e1baa548cf01))
- **visualizations:** add view async capapbility ([d2eb414](https://github.com/biosimulations/Biosimulations/commit/d2eb414a45e9d4eee2caf0edcf87e29c1b066add))

### Reverts

- "Formatted Files. [skip ci]" ([0414aca](https://github.com/biosimulations/Biosimulations/commit/0414aca53ff1c5ee5de9920fdb3e7e5810582a1b))
- Revert "debug redis host" ([7264f87](https://github.com/biosimulations/Biosimulations/commit/7264f87a60c994a655fd80536c33eb234efb5d2b))
- Revert "Feat(Combine-serive): Update API Specification" ([7ada4ce](https://github.com/biosimulations/Biosimulations/commit/7ada4ce33538e4bcd3bae9ee798ea46e417cee99))
- Revert "Updated Ontologies" ([aab6f70](https://github.com/biosimulations/Biosimulations/commit/aab6f708ab3f2d6e39f780bcb4bff1a188376748))
- Revert "Bump @sendgrid/mail from 7.4.1 to 7.4.2 in /biosimulations (#1975)" (#1979) ([a96e0ea](https://github.com/biosimulations/Biosimulations/commit/a96e0ea8f68b97fc2722d464614ee002bcd479c5)), closes [#1975](https://github.com/biosimulations/Biosimulations/issues/1975) [#1979](https://github.com/biosimulations/Biosimulations/issues/1979)
- Revert "Bump @nrwl/cli from 10.4.1 to 11.0.20 in /biosimulations (#1855)" (#1860) ([b7d637e](https://github.com/biosimulations/Biosimulations/commit/b7d637ec1decf056c1642608d821c7b53422b46d)), closes [#1855](https://github.com/biosimulations/Biosimulations/issues/1855) [#1860](https://github.com/biosimulations/Biosimulations/issues/1860)
- Revert "merging enumerations of simulation status; aligning names of properties 'resultsSize' and 'resultSize'" ([08e96f2](https://github.com/biosimulations/Biosimulations/commit/08e96f2fe8c7e88c7fbb2a6653416182fc438700))
- Revert "Revert "styling nested lists and lists after paragraphs"" ([16a2d5b](https://github.com/biosimulations/Biosimulations/commit/16a2d5ba5e39a79c64bca5af77c1382c76dfb2c5))
- Revert "Revert "adding management for app-specific configuration (e.g,, appName, logo, etc."" ([1c6190a](https://github.com/biosimulations/Biosimulations/commit/1c6190ab50c95359233e04befc2c778fa1dbd5c1))
- Revert "Revert "adding documentation"" ([f71cfcd](https://github.com/biosimulations/Biosimulations/commit/f71cfcddcf15006f3f84660517c7fdf52ada0dbe))
- Revert "Revert "editing help"" ([4aea612](https://github.com/biosimulations/Biosimulations/commit/4aea6128b0f4fbb2d881efbcc5eb2416527f05c9))
- Revert "Revert "adding documentation of supported SED-ML features"" ([32beb34](https://github.com/biosimulations/Biosimulations/commit/32beb34d5a823e3975d6115d8f225044fe6297a5))
- Revert "changing completed to updated" ([d9fe18e](https://github.com/biosimulations/Biosimulations/commit/d9fe18edbbf4a1d0307caecdeb7a6927b0338abf))
- Revert "update to angular 10" ([e1b9fcf](https://github.com/biosimulations/Biosimulations/commit/e1b9fcfae21f5bd8b19e4cee23a8df2b82c1df02))
- Revert "Bump husky from 4.0.1 to 4.0.5 in /CRBM-Viz (#290)" ([4dd2cc3](https://github.com/biosimulations/Biosimulations/commit/4dd2cc3eee2fcf33bb24ca8ef999162823e875a7)), closes [#290](https://github.com/biosimulations/Biosimulations/issues/290)

- feat (dispatch-api): remove download endpoint ([785ad27](https://github.com/biosimulations/Biosimulations/commit/785ad27f1477ac6122c2a735cb46201928a0f754))
- feat (dispatch) : Change datamodel of returned results ([ba42dcc](https://github.com/biosimulations/Biosimulations/commit/ba42dcc8b74068a280c6dc2f7d915c8c27a55f45))

### BREAKING CHANGES

- The /download endpoint has been removed. Should be replaced by /results/download
- The results are now returned as an array of objects (AOS) rather than an object of arrays (SOA)
# Source code

The source code for the BioSimulations and BioSimulators platforms is available from [our Git repository](https://github.com/biosimulations/biosimulations). The source code for the underlying simulation capabilities is available from the [BioSimulators GitHub organization](https://github.com/biosimulators).
# Guide to contributing to the BioSimulations and BioSimulators platforms

We enthusiastically welcome contributions to BioSimulations and BioSimulators! This document describes how developers can contribute the BioSimulations and BioSimulators platforms (e.g., web applications, REST APIs, databases, and simulation services). Information for investigators about contributing simulation projects to BioSimulations is available [here](../users/publishing-projects.md). Information for simulation software tool developers about contributing simulation tools to BioSimulators is available at [here](../users/publishing-tools.md).

## Coordinating contributions

Before getting started, please contact the lead developers at [info@biosimulations.org](mailto:info@biosimulations.org) to coordinate your planned contributions with other ongoing efforts. Please also use GitHub issues to announce your plans to the community so that other developers can provide input into your plans and coordinate their own work. As the development community grows, we will institute additional infrastructure as needed, such as a leadership committee and regular online meetings.

## Repository organization

The repository is organized as a monorepo using tooling from [Nx](https://nx.dev/angular/getting-started/why-nx). A brief guide to working with Nx is available [here](./setup/nx-tutorial.md).

### Apps

The `apps` directory contains code for the high level applications that are a part of BioSimulations and BioSimulators. An application is code that is deployed independently, either as a website, API, or backend service. The applications are implemented using TypeScript, except the COMBINE API, which is implemented using Python. See the [README](https://github.com/biosimulations/biosimulations/blob/dev/apps/combine-api/README.md) for the COMBINE API for more information about installing, running, and editing the COMBINE API.

### Libraries

The `libs` directory contains code that can be used by multiple applications. Each library contains a `README.md` file that describes its purpose. The libraries are implemented using TypeScript. To enforce proper separation of concerns and manage dependency trees, BioSimulations and BioSimulators employ constraints on the libraries that can be used by each application. For more information, read the [Nx documentation](https://nx.dev/angular/workspace/structure/monorepo-tags), and look at the repository's [linting rules](https://github.com/biosimulations/biosimulations/blob/dev/.eslintrc.json). In general, BioSimulations and BioSimulators follow the [principles recommended by Nx](https://nx.dev/angular/guides/monorepo-nx-enterprise).

## Coding convention

BioSimulations and BioSimulators follow standard TypeScript style conventions:

- Module names: `lower-dash-case`
- Class names: `UpperCamelCase`
- Function names: `lowerCamelCase`
- Variable names: `lowerCamelCase`

Further information about style conventions can be found in the [lint rules definition](https://github.com/biosimulations/biosimulations/blob/dev/.eslintrc.json).

## Linting

We strive to create high-quality code. BioSimulations and BioSimulators use [ESLint](https://eslint.org/) to identify potential errors. ESLint can be executed by running `nx affected:lint` to lint all applications and libraries or `nx run {app-or-lib}:lint` to lint an individual application or library.

See the [README](https://github.com/biosimulations/biosimulations/blob/dev/apps/combine-api/README.md) for the COMBINE API for more information about linting the COMBINE API.

## Testing

We strive to have complete test coverage of BioSimulations and BioSimulators. As such, all contributions to BioSimulations and BioSimulators should be tested. In particular, each module should be accompanied by a specification test in the same directory with the extension `.spec.ts`. All applications and libraries can be tested by running `nx affected:test`. Individual applications or libraries can be tested by running `nx run {app-or-lib}:test`.

Upon each push and pull request to GitHub, GitHub will trigger actions to execute all of the tests.

See the [README](https://github.com/biosimulations/biosimulations/blob/dev/apps/combine-api/README.md) for the COMBINE API for more information about testing the COMBINE API.

## Submitting changes

Please use GitHub pull requests to submit changes. Each request should include a brief description of the new and/or modified features. Upon each pull request, GitHub will trigger actions to lint and test all of the applications and libraries. Pull requests will be approved once all tests are passing and the pull request is reviewed by one of the lead developers.

### Commit convention

BioSimulations and BioSimulators uses [conventional commits](https://www.conventionalcommits.org/) to enable changelog generation and versioning. Developers are encouraged to use [Commitizen](http://commitizen.github.io/cz-cli/) to easily create compliant commit messages. To use Commitizen, simply run `npm run commit` instead of ` git commit` to commit changes.

## Releasing and deploying new versions

Contact [info@biosimulations.org](mailto:info@biosimulations.org) to request release and deployment of changes.
# Code of conduct for developers of the BioSimulations and BioSimulators platforms

## Our pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our standards

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

## Our responsibilities

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
reported by contacting the lead developers at [info@biosimulations.org](mailto:info@biosimulations.org). All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [https://www.contributor-covenant.org/version/1/4/code-of-conduct.html](https://www.contributor-covenant.org/version/1/4/code-of-conduct.html).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
[https://www.contributor-covenant.org/faq](https://www.contributor-covenant.org/faq).
# Architectural decisions and philosophy

This page outlines some of the rationale for the architecture of the project.

## Hybrid cloud-on-premise architecture

BioSimulations uses a hybrid architecture in which the APIs, backend services, and databases are deployed in a Kubernetes cluster in the commercial cloud; simulation runs are executed using on-premise high-performance computing (HPC) at UConn Health; and simulation projects and their results are stored using on-premise resources at UConn Health. We chose this hybrid architecture to minimize the long-term cost of deploying BioSimulations, particularly executing simulations for the community. Nevertheless, BioSimulations and BioSimulators are architected such they can be deployed to alternative resources with minimal modification.
# Deploying BioSimulations in alternative environments

As discussed [here](./philosophy.md), although BioSimulations use Slurm to run simulations, BioSimulations is designed to be easily modified to execute simulations with other environments.

Simulations are executed by the [`dispatch-service`](https://github.com/biosimulations/biosimulations/tree/dev/apps/dispatch-service/src). This service generates Slurm scripts to run simulations, submits these scripts via an SSH connection to a Slurm queue, and uses an SSH connection to monitor simulation jobs. When jobs complete, the service retrieves their logs via SSH and submits them to BioSimulations' API. 

The following modifications are needed to deploy BioSimulations in other systems:

- Slurm-based HPCs: The path to the `dispatch-service` user's home directory, the location of the simulation file, and other configuration variables must be modified. This is done simply by setting the appropriate environment variables through the application config. An example configuration is available [here](https://github.com/biosimulations/deployment/blob/main/config/prod/dispatch-service/config.env).

- Other job schedulers: The generation of the Sbatch script and submission of the job must be modified to support the desired scheduler. The job submission and monitoring code must also be modified to support the desired scheduler. The following files contain the code that must be modified:
    - [`sbatch.service.ts`](https://github.com/biosimulations/biosimulations/blob/dev/apps/dispatch-service/src/app/services/sbatch/sbatch.service.ts)
    - [`hpc.service.ts`](https://github.com/biosimulations/biosimulations/blob/dev/apps/dispatch-service/src/app/services/hpc/hpc.service.ts)

- Cloud-based compute: Using a cloud-based environment may also require changing how simulation runs are monitored and how their logs are retrieved. The following files contain the code that may need to be modified:
    - [`monitor.processor.ts`](https://github.com/biosimulations/biosimulations/blob/dev/apps/dispatch-service/src/app/submission/monitor.processor.ts)
    - [`complete.processor.ts`](https://github.com/biosimulations/biosimulations/blob/dev/apps/dispatch-service/src/app/submission/complete.processor.ts)

For further assistance, please contact the [BioSimulations Team](mailto:info@biosimulations.org).
# Technological foundation
BioSimulations and BioSimulators are implemented using several open-source tools and cloud platforms, including those highlighted below.

<div class="logos">
<div class="logos-row">
    <a href="https://angular.io" rel="noopener" target="_blank" title="Angular">
    <img class="zoom" src="/assets/images/about/tech/angular.svg" />
    </a>

    <a href="https://auth0.com/" rel="noopener" target="_blank" title="Auth0">
    <img class="zoom" src="/assets/images/about/tech/auth0.svg" />
    </a>

    <a href="https://www.docker.com/" rel="noopener" target="_blank" title="Docker">
    <img class="zoom" src="/assets/images/about/tech/docker.svg" />
    </a>

    <a
    href="https://realfavicongenerator.net/"
    rel="noopener" target="_blank"
    title="Favicon Generator"
    >
    <img
        class="zoom"
        src="/assets/images/about/tech/real-favicon-generator.png"
    />
    </a>

    <a
    href="https://fontawesome.com"
    rel="noopener" target="_blank"
    title="Font Awesome"
    >
    <img class="zoom" src="/assets/images/about/tech/font-awesome.svg" />
    </a>

    <a href="https://www.github.com/" rel="noopener" target="_blank" title="GitHub">
    <img class="zoom" src="/assets/images/about/tech/github.svg" />
    </a>

    <a href="https://gravatar.com/" rel="noopener" target="_blank" title="Gravatar">
    <img class="zoom" src="/assets/images/about/tech/gravatar.svg" />
    </a>

    <a href="https://kubernetes.io/" rel="noopener" target="_blank" title="Kubernetes">
    <img class="zoom" src="/assets/images/about/tech/kubernetes.svg" />
    </a>

    <a href="https://www.mongodb.com/" rel="noopener" target="_blank" title="MongoDB">
    <img class="zoom" src="/assets/images/about/tech/mongodb.svg" />
    </a>

    <a
    href="https://material.io/"
    rel="noopener" target="_blank"
    title="Material Design"
    >
    <img class="zoom" src="/assets/images/about/tech/material.svg" />
    </a>

    <a
    href="https://materialtheme.arcsine.dev/"
    rel="noopener" target="_blank"
    title="Material Theme Generator"
    >
    <img
        class="zoom"
        src="/assets/images/about/tech/material-theme-generator.svg"
    />
    </a>

    <a rel="noopener" target="_blank" href="https://nestjs.com/" title="NestJS">
    <img class="zoom" src="/assets/images/about/tech/nestjs.svg"
    /></a>
    <a rel="noopener" target="_blank" href="https://www.netlify.com" title="Netlify">
    <img class="zoom" src="/assets/images/about/tech/netlify.svg"
    /></a>
    <a rel="noopener" target="_blank" href="https://www.openapis.org/" title="OpenAPI">
    <img class="zoom" src="/assets/images/about/tech/openapi.svg"
    /></a>
    <a rel="noopener" target="_blank" href="https://spdx.org" title="SPDX">
    <img class="zoom" src="/assets/images/about/tech/spdx.svg"
    /></a>
    <a rel="noopener" target="_blank" href="https://swagger.io" title="Swagger">
    <img class="zoom" src="/assets/images/about/tech/swagger.svg"
    /></a>
</div>
</div>
# Architecture

The BioSimulations/BioSimulators platform is a distributed computing system that is architected for extensibility, scalability, and ease of development. By separating the various components of the platform, we hope to enable developers to focus on relevant portions of the system and not on the details of the underlying architecture. In deciding which parts of the system should be separated, we try to balance the added complexity of distributed computing with clean separation of concerns. Details about the decisions and thought processes behind the architecture can be found in [Architecture philosophy](./philosophy.md) page.

The components of the platform are organized as follows:

- Front-end applications
    - [BioSimulations](https://biosimulations.org) (`platform` app)
    - [runBiosimulations](https://run.biosimulations.org) (`dispatch` app)
    - [BioSimulators](https://simulators.org) (`simulators` app)
- Public APIs
    - [BioSimulations API](https://api.biosimulations.org) (`api` app)
    - [BioSimulators API](https://api.biosimulators.org) (`simulators-api` app)
- Back-end services
    - Dispatch service (`dispatch-service` app)
    - [COMBINE API](https://combine.api.biosimulations.org) (`combine-api` app)
- Data storage
    - Mongo database
    - [S3 buckets](https://files.biosimulations.org) (accessed through `api` app)
- Computing infrastructure
    - GKE Kubernetes cluster
    - High performance computing cluster managed with SLURM
# Running required backend services for local development

## Background Services

!!! tip 
    We recommend using [VsCode](https://code.visualstudio.com/) for developing BioSimulations. The BioSimulations git repository contains a [development container](https://code.visualstudio.com/docs/remote/containers) [configuration](https://github.com/biosimulations/biosimulations/blob/dev/.devcontainer/devcontainer.json) that simplifies the setup of the local environment.

The BioSimulations apps requires connecting to several infrastructure services for functions such as messaging, database, and storage. Developing BioSimulations locally requires access to the following services:

- Redis

    [Redis](https://redis.io/) is used for caching the results of the API and for managing queues for submitting, monitoring, and processing simulation runs.

    We recommend running a local redis container with the following command:

    ```bash
    docker run -d -p 6379:6379 --network host --name redis redis
    ```

- NATS messaging queue

    [NATS](https://docs.nats.io/) is used for messaging between the API and the backend services. Currently, the NATS connection does not support guaranteed message delivery semantics. The latest [NATS Jetstream](https://docs.nats.io/nats-concepts/jetstream) adds support for exactly-once delivery semantics, which will be can be for submission jobs in the future.

    We recommend running a local NATS container with the following command:

    ```bash
    docker run -d -p 4222:4222 --network host --name nats nats
    ```

- MongoDB

    [MongoDB](https://docs.mongodb.com/) is used as our primary database and contains information about the simulation runs, their logs, specifications, metadata, etc. It is accessed through the API.

    We recommend running a local MongoDB container with the following command:
    ```bash
    docker run -d -p 27017:27017 --network host --name mongodb mongo
    ```

    Alternatively, you can use a free MongoDB cluster from [MongoDB Atlas](https://www.mongodb.com/cloud/atlas/).

    These methods will allow you to develop BioSimulations locally, but will not have access to simulation stored on the development server. If you would like to develop against the dev biosimulations database, please [contact us](/about/contact/) for access.

- S3

    We recommend using an [AWS S3](https://aws.amazon.com/s3/) or [Google Cloud Storage](https://cloud.google.com/storage/) bucket for local development. If you would like to develop against the dev biosimulations bucket, please [contact us](/about/contact/) for access.

    The dev biosimulations bucket can be accessed in the following ways: 
    
    - Connect to the UCONN VPN and use the endpoint `https://s3low.scality.uchc.edu`

    - Connect to the BioSimulations submit node on the UCONN HPC via SSH and use local port bindings. The bucket can then be accessed via 'https://localhost:4443.' You will likely need to disable ssl validation in the library, tool, or code that you are using to connect to the bucket. 

    ```bash
     ssh -i ~/.ssh/id_hpc -L localhost:4443:s3low.scality.uchc.edu:443 crbmapi@biosim-submit-ext.cam.uchc.edu
    ```
     Where `id_hpc` is a private key that has access to the biosimulations user account. 

     To request a ssh key, please [contact us](/about/contact/).

- HPC
    
    Accessing the an appropriate HPC backend is the most bespoke part of setting up a local environment for development on the BioSimulations platform. If your development work does not require running simulations on the HPC, we recommend temporarily modifying the code to skip submitting jobs to the HPC or returning mock responses. This will allow you to continue to develop the API, and other aspects of the system without access to the HPC. For more information about which parts of the system interact with the HPC, see the [architecture deployment information](../architecture/deployment.md).

     If you specifically require access to the HPC, please [contact us](/about/contact/).


# Nx tutorial

This project uses [Nx](https://nx.dev). Below is a brief introduction to using Nx with this project.

## Generate a skeleton for an application

Run `ng g @nrwl/angular:app my-app` to generate an application.

When using Nx, you can create multiple applications and libraries in the same workspace.

## Generate a skeleton for a library

Run `ng g @nrwl/angular:lib my-lib` to generate a library.

Libraries are sharable across libraries and applications. They can be imported as `@biosimulations/mylib`.

## Generate a skeleton for a component

Run `ng g component my-component --project=my-app` to generate a skeleton for a new component.

## Run a development server

Run `ng serve my-app` to run a development server. Then navigate to [http://localhost:4200/](http://localhost:4200/). The application will automatically reload if you change any of the source files.

## Run unit tests for an application or library

Run `ng test my-app` to execute the unit tests via [Jest](https://jestjs.io).

Run `nx affected:test` to execute the unit tests affected by your recent changes.

## Run end-to-end tests for an application or library

Run `ng e2e my-app` to execute the end-to-end tests via [Cypress](https://www.cypress.io).

Run `nx affected:e2e` to execute the end-to-end tests affected by your recent changes.

## Build an application

Run `ng build my-app` to build the project. The build artifacts will be stored in the `dist/` directory. Use the `--prod` flag for a production build.

## Mapping the dependencies of the applications and libraries

Run `nx dep-graph` to see a diagram of the dependencies of your projects.

## Cloud-based computation memoization

This project uses [Nx Cloud](https://nx.app/) to help developers save time when building and testing. Visit [Nx Cloud](https://nx.app/) to learn more.

## Further help

Below are several resources for additional help:

* [10-minute video which outlines the features of Nx](https://nx.dev/angular/getting-started/what-is-nx)
* [Interactive tutorial](https://nx.dev/angular/tutorial/01-create-application)
* [Documentation](https://nx.dev/angular)
# Configuring endpoints for local environments
## Loading endpoints used by the applications
The endpoints for the front end and backend applications are located in the shared configuration library, located at `libs/config/shared`. The endpoints are loaded dynamically depending on the value of the `env` parameter provided when initializing the `Endpoints` class. 

The loading of the endpoints differs slightly depending on whether the application is running a frontend (browser) or backend (server) application.
### Frontend Applications
For front end applications, the default value of the `env` parameter is read from the `@biosimulations/config/shared` library, located at `libs/shared/environments`. This library contains several different environment definitions such as "development", "staging", "production", etc. The `environment.ts` file loads one of these definitions.  The `env` parameter that is loaded is then fed into the `Endpoints` class as described above.

To configure which endpoints are loaded, change the name of the file being loaded in the `environment.ts` file. For example, to load the endpoints for the "local" environment, change the name of the file to `environment.local`.

```typescript
import { environmentType } from './environment.type';
// Change the name of the file to environment.type.ts where type is the name of the environment you wish to load
import { environment as currentEnvironment } from './environment.dev';
export const environment: environmentType = currentEnvironment;
```

### Backend Applications

For back end applications, the developer must provide the value of the `env` parameter when initializing the `Endpoints` class. This value is then used to load the appropriate endpoints.

In most cases, the value of the `env` parameter should be loaded from the the configuration service provided by the `@biosimulations/config/nest` library located at `libs/config/nest`. The following is an example of how to load the endpoints for the current environment.

```typescript
import { Injectable } from '@nestjs/common';
import { Endpoints } from '@biosimulations/config/common';
import { ConfigService } from '@nestjs/config';

@Injectable()
export class SomeService {
  private endpoints: Endpoints;

  public constructor(
    private configService: ConfigService,
  ) {
    const env = configService.get('server.env');
    this.endpoints = new Endpoints(env);
  }
}
```
Unlike, the front end applications, the backend applications are capable of loading endpoints dynamically. The `Endpoints` class loads overrides the default endpoints for the current environment by loading environment variables. In order to override a specific endpoint, set the environment variable with the name of the endpoint to the value you wish to use. For example, to override the `combine-api` endpoint, set the environment variable `COMBINE_API_URL` to the value you wish to use.
A list of the environment variables that can be overridden is located in the `EndpointLoader` class of the `@biosimulations/config/shared` library.
### External Endpoints 

Due to the distributed architecture of the application, various endpoints may not be accessible from different locations. For example, if a developer is developing their application locally, they may set the `API_URL` environment to `http://localhost:3333` to allow the locally running `dispatch-service` to post data to their locally running API. However, the HPC cluster would not be able to download the simulation project to execute from a `localhost` address, since this is not accessible from the HPC cluster. 

For this reason, the `Endpoints` class provides an `external` variant of each endpoint. This is the endpoint that should be used when urls are being shared outside of the 'current local' environment. The application developers must be sure to load these endpoints from the `Endpoint` class at the appropriate points in the code. These endpoints can be overridden the same way as the default endpoints, but adding the `EXTERNAL` prefix. For example, to to address the above issue, the user would set the `EXTERNAL_API_URL` environment variable to an address that is accessible from the HPC cluster, such as `https://api.biosimulations.dev`, or a url provided by an tunneling service such as [localtunnel](http://localtunnel.github.io/www/).

!!! warning 
    Make sure you understand the security implications of exposing locally running applications on your machine to the world via public urls. In particular, make sure that these URLs are not accidentally committed to the repository.# Security policy

## Introduction

BioSimulations and BioSimulators welcome feedback from security researchers and the general public to help improve our security. If you believe you have discovered a vulnerability, privacy issue, exposed data, or other security issues in any of our assets, we want to hear from you. This policy outlines steps for reporting vulnerabilities to us, what we expect, and what you can expect from us.

## Systems in scope

This policy applies to any assets owned, operated, or maintained by BioSimulations or BioSimulators.

BioSimulations and BioSimulators are updated on a rolling release. If any security issues are found, only the latest releases and deployments will contain any security fixes.

## Systems out of scope

This policy does not apply to assets not owned, operated, or maintained by BioSimulations or BioSimulators. 

Vulnerabilities discovered or suspected in out-of-scope systems should be reported to the appropriate vendor or applicable authority.

## Our commitments

When working with us, you can expect us to:

- Respond to your report promptly, and work with you to understand and validate your report;
- Work to remediate discovered vulnerabilities in a timely manner, within our operational constraints;
- Strive to keep you informed about our progress addressing the vulnerability; and
- Extend Safe Harbor for your vulnerability research that is related to this policy.

## Our expectations

In participating in our vulnerability disclosure program in good faith, we ask that you:

- Play by the rules, including following this policy and any other relevant agreements. If there is any inconsistency between this policy and any other applicable terms, the terms of this policy will prevail;
- Report any vulnerability you’ve discovered promptly;
- Avoid violating the privacy of others, disrupting our systems, destroying data, and/or harming user experience;
- Use only our development deployments (eg. api.biosimulations.dev) for testing and research purposes unless you strongly believe there are deployment-specific vulnerabilities;
- Use only the official channels outlined below to discuss vulnerability information with us;
- Provide us a reasonable amount of time (at least 90 days from the initial report) to resolve the issue before you disclose it publicly;
- Perform testing only on in-scope systems, and respect systems and activities which are out-of-scope;
- If a vulnerability provides unintended access to data: Limit the amount of data you access to the minimum required for effectively demonstrating a proof of concept; and cease testing and submit a report immediately if you encounter any user data during testing, such as Personally Identifiable Information (PII), Personal Healthcare Information (PHI), credit card data, or proprietary information;
- You should only interact with test accounts you own or with explicit permission from the account holder; and
- Do not engage in extortion.  

## Official channels for reporting vulnerabilities

Please report security issues to our lead developers via [info@biosimulations.org](mailto:info@biosimulations.org). Please provide as much relevant information as possible. The more details you provide, the easier it will be for us to triage and fix the issue.

## Safe harbor

When conducting vulnerability research, according to this policy, we consider this research conducted under this policy to be:

- Authorized concerning any applicable anti-hacking laws, and we will not initiate or support legal action against you for accidental, good-faith violations of this policy;
- Authorized concerning any relevant anti-circumvention laws, and we will not bring a claim against you for circumvention of technology controls;
- Exempt from restrictions in our Terms of Service (TOS) that would interfere with conducting security research, and we waive those restrictions on a limited basis; and
- Lawful, helpful to the overall security of the internet, and conducted in good faith.

You are expected, as always, to comply with all applicable laws. If legal action is initiated by a third party against you and you have complied with this policy, we will take steps to make it known that your actions were conducted in compliance with this policy.

If at any time you have concerns or are uncertain whether your security research is consistent with this policy, please submit a report through our official channel before going any further.

!!!note
    Safe Harbor applies only to legal claims under the control of BioSimulations and BioSimulators. This policy does not bind independent third parties.
# Contact Us

## :thought_balloon: Questions, comments, or feedback 

We welcome any comments, questions, or discussion about the project. Please create a discussion or question in our [discussion forum](https://github.com/biosimulations/biosimulations/discussions).

## :envelope: Contact the developers individually 
To privately contact the BioSimulations team, you can send us an email at [info@biosimulations.org](mailto:info@biosimulations.org).
# Core development team

BioSimulations and BioSimulators were developed by the [Center for Reproducible Biomedical Modeling](http://reproduciblebiomodels.org/) including [Bilal Shaikh](https://bshaikh.com) and [Jonathan Karr](https://www.karrlab.org) at the [Icahn School of Medicine at Mount Sinai](https://icahn.mssm.edu); Akhil Marupilla, [Dan Vasilescu](https://www.linkedin.com/in/dvasilescu), [Mike Wilson](https://www.linkedin.com/in/mike-wilson-08b3324/), [Michael Blinov](https://health.uconn.edu/blinov-lab/), and [Ion Moraru](https://facultydirectory.uchc.edu/profile?profileId=Moraru-Ion") at the [Center for Cell Analysis and Modeling](https://health.uconn.edu/cell-analysis-modeling/) at UConn Health; Lucian Smith and [Herbert Sauro](https://www.sys-bio.org/) at the University of Washington.

<div class="logos">
    <div class="logos-row">
        <a
        href="https://reproduciblebiomodels.org/"
        rel="noopener" target="_blank"
        title="Center for Reproducible Biomedical Modeling"
        >
            <img class="zoom" src="/assets/images/about/team/crbm.svg" />
        </a>

        <a
            href="https://www.karrlab.org/"
            rel="noopener" target="_blank"
            title="Center for Reproducible Biomedical Modeling"
        >
            <img class="zoom" src="/assets/images/about/team/karr-lab.svg" />
        </a>

        <a
            href="https://icahn.mssm.edu"
            rel="noopener" target="_blank"
            title="Icahn School of Medicine at Mount Sinai"
        >
            <img class="zoom" src="/assets/images/about/team/sinai.svg" />
        </a>

        <a
            href="https://health.uconn.edu/"
            rel="noopener" target="_blank"
            title="UConn Health"
        >
            <img class="zoom" src="/assets/images/about/team/uconn.svg" />
        </a>

        <a
            href="https://uw.edu"
            rel="noopener" target="_blank"
            title="University of Washington"
        >
            <img class="zoom" src="/assets/images/about/team/uw.svg" />
        </a>
    </div>
</div>
# Terms of service

## Permitted usage
BioSimulations, runBioSimulations, and BioSimulators promote open science by providing the community a central platform for sharing and executing biomodeling projects and simulation tools. The public simulation projects, simulation tools, and simulation results available through BioSimulations and BioSimulators are available under the open-source licenses chosen by their contributors. These licenses are noted with each resource. BioSimulations imposes no additional restrictions on the use of these community-contributed resources beyond those imposed by their contributors. These licences describe who can use each resource, for which purposes, and what attribution may be required.

To publish a simulation tool, execute a simulation project, or publish a simulation project, users must agree to grant BioSimulators, runBioSimulations, and BioSimulators free of charge, a license to store the tool, project, and/or simulation results.

## Attribution
BioSimulations and BioSimulators request attribution (e.g., in publications, services, or products) for use of their simulation projects, simulation tools, and online services in accordance with good scientific practice. Please use the citations indicated [here](./citations.md).

## Privacy policy
BioSimulations and BioSimulators collect limited personal data. Our [privacy policy](./privacy.md) describes how we collect, store, and use this data and who has access to this data.

## Disclaimers
The simulation tools available through runBioSimulations and BioSimulators and simulation results provided by runBioSimulations and BioSimulations are provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the BioSimulations or BioSimulators Teams be liable for any claim, damages or other liability whatsoever, whether in an action of contract, tort or otherwise, arising from, out of or in connection with these resources or the use or other dealings in these resources.

The BioSimulations and BioSimulators Teams do not guarantee the accuracy of the simulation tools available through BioSimulators or runBioSimulations or the simulation results provided by runBioSimulations and BioSimulations, nor the suitability of these resources for any purpose.

The BioSimulations and BioSimulators Teams will make reasonable effort to maintain continuity of BioSimulations, runBioSimulations, and BioSimulators. However, the BioSimulations and BioSimulators Teams accept no responsibility for the consequences of any temporary or permanent discontinuity in service.

Any feedback provided to the BioSimulations and BioSimulators Teams will be treated as non-confidential unless the individual or organization providing the feedback states otherwise.

Any attempt to use BioSimulations or BioSimulators to a level that prevents, or looks likely to prevent, BioSimulations or BioSimulators from servicing other users may result in the use being blocked.

If you post offensive, inappropriate or objectionable content to BioSimulations or BioSimulators, or otherwise engage in disruptive behavior, we may use information from our user behavior logs to stop such behavior. Where we believe that you are, or may be in breach, of any applicable laws, we may use information from our user behavior logs to inform relevant third parties about the content and your behavior.

The BioSimulations and BioSimulators Teams do not accept responsibility for personal data voluntarily disclosed through files submitted to BioSimulations or BioSimulators, nor for the consequences of any breach of the confidentiality of simulation studies stored in runBioSimulations or BioSimulations by third parties.

The BioSimulations and BioSimulators Teams are not liable to you or third parties claiming through you, for any loss or damage.

## Updates to this policy
We reserve the right to update this policy at any time. The date of the most recent revision of this policy is at the bottom of this page.

## Contact
If you have any questions about this policy, please contact Jonathan Karr and Ion Moraru at [info@biosimulations.org](mailto:info@biosimulations.org) or at 263 Farmington Avenue, Farmington, CT 06030-6406, USA.
# Privacy policy

## Personal data that we collect
Optionally, we collect profiles of the owners of projects published to BioSimulations and simulation tools published to BioSimulators. This includes their name, institution, website, email, social profiles (e.g., ORCID), brief biography, and head shot. Optionally, we also collect email addresses from users of runBioSimulations. In addition, we collect standard information about user behavior, such as the operating system and browser used to access our websites and the dates and times that users visited our websites.

## How we use personal data
We use the author profiles to publicly credit authors for their projects and simulation tools. We use the email addresses to notify users when their simulations have completed and are ready for analysis. We use the user behavior data to help improve the functionality of BioSimulations and BioSimulators for the scientific community. This includes helping us better understand the needs of users.

## How we collect personal data
After logging in, authors will be able to edit their profile using our online form. Users can optionally provide their email addresses when they submit simulations through the simulation submission form or API endpoint. We use performance cookies and Google Analytics to collect data about user behavior. We do not use tracking cookies.

## How we store personal data
The author profiles and user email addresses are stored in our MongoDB database hosted by MongoDB Atlas. We use Google Analytics to store the user behavior data. Personal data will not be transferred to international organizations.

## How long we keep personal data
This personal data will be stored as long as BioSimulations and BioSimulators are live.

## Who has access to this personal data
The author profiles are publicly accessible. The email addresses and user behavior data are only accessible to the lead developers. We may also make information about the total volume of usage of particular simulation projects and simulation tools available to the public and third party users and organizations who supply the projects and tools without details of any individual's use.

## Your rights regarding your personal data
You have the right to view, edit, or delete your author profile at any time. Because we do not store identifying user behavior data, we cannot provide information about the user behavior data processed about you.

## How to contact us
If you have any questions about this policy, please contact Jonathan Karr and Ion Moraru at [info@biosimulations.org](mailto:info@biosimulations.org) or at 263 Farmington Avenue, Farmington, CT 06030-6406, USA.
# Our use of cookies and other similar local storage

## What are cookies and other local storage?

Cookies and local storage are mechanisms for storing information related to your use of BioSimulators and BioSimulations in your browser.

More information about cookies and local storage is available from [allaboutcookies.org](https://allaboutcookies.org) and [gdpr.eu/cookies](https://gdpr.eu/cookies).

## How do we use cookies and local storage?

BioSimulators and BioSimulations use cookies and local storage to improve your experience on our websites. This includes:

- Providing essential functionality, such as authentication and security.
- Helping us understand how users use our websites. This information helps us improve our websites.
- Collecting anonymous, aggregated information about visitors to our websites. This information helps us plan the resources needed to run our websites and report their impact to the funding agencies that support BioSimulators and BioSimulations.
- Enabling users to save their settings, such as their cookie preferences.

## What types of cookies and local storage do we use?

We only use first-party cookies and local storage. This means that we do not share cookies with any third parties or external companies. 

Our cookies and local storage do not contain any personal information. However, we do generate cookies that are unique to each client. These cookies cannot be used to identify you.

Below is a list of the types of cookies and local storage that we use. More information about each type of cookie is available at [allaboutcookies.org](https://allaboutcookies.org) and [gdpr.eu/cookies](https://gdpr.eu/cookies).


### Necessary cookies and local storage
These cookies are essential for our websites to function. 

Currently, we use cookies and local storage for the following purposes:

- Store cookie preferences: When you visit our websites and select your cookie preferences, we store those selections in your browser's local storage. This ensures that we do not collect information that you have not explicitly allowed.

- Enable authentication: For sections of our websites that require you to log into an account, we use cookies to verify that you are logged in. This enables us to ensure that only you have access to your account and its private resources. We use the following cookies for this purpose:
    - `auth0`, `auth0_compat`: These cookies are used to verify that you are logged in.
    - `did`, `did-compat`: These cookies are used to protect BioSimulators and BioSimulations from bots and other automated attacks.

- Store information about the simulation runs that you have submitted in the browser's local storage, to enable you to view the simulation runs you have loaded and submitted.

### Functional cookies and local storage

These cookies are used to support additional, non-essential functionality that provides users with a more personalized experience.

- Currently, we do not use functional cookies


### Performance cookies and local storage

These cookies are used to collect data for improving the performance of our websites.

Currently, we use cookies for the following purposes:

- `biosim_analytics_ga_XXXXXXXXXX`: This cookie allows us to distinguish between users who are new to our websites and users who have visited our websites before. This helps us accurately measure the utilization and performance of our websites.

- `biosim_analytics_la_XXXXXXXXXX`: This cookie is used to monitor how users interact with our websites. For example, which links users follow to navigate our websites. This information helps us improve our websites for users.

These cookies are used to gather information that is sent to Google Analytics. No personal information is collected or sent to Google Analytics. For more information, see <!-- * #no-spell-check* --> [Google Analytics](https://developers.google.com/analytics/devguides/collection/gtagjs/cookie-usage). 

### Tracking cookies and local storage

These cookies are used to track user's behavior and history for marketing and advertising purposes.

BioSimulations and BioSimulators DO NOT USE any tracking or targeting cookies for any purpose. 

BioSimulations and BioSimulators DO NOT COLLECT any personal information through the use of cookies.

## How to manage your cookies and local storage

You can set your browser to accept cookies, to see what cookies and local storage are being set for you, or to delete cookies or local storage.

You can also set your browser to block cookies. However, this may affect the functionality of our websites.

Upon your first visit to our websites, we will present you with an explanation of the cookies we use and options for disabling some or all of them. These selections will be saved in your browser's local storage.

If you wish to change your cookie preferences, clear the site data in your browser for our websites. When you visit our websites again, you will again be prompted to select your cookie preferences. More information about clearing your site data is available at [allaboutcookies.org](https://allaboutcookies.org) and [gdpr.eu/cookies](https://gdpr.eu/cookies).

## Updates to this policy

As we continue to develop BioSimulations and BioSimulators, we may have additional needs or technical requirements that may require us to change our this policy. We will post any changes to this policy on this page. Below is the date when this policy was last updated.
# Funding

BioSimulations and BioSimulators were developed with support from the [Center for Reproducible Biomodeling Modeling](https://reproduciblebiomodels.org) from the [National Institute of Bioimaging and Bioengineering](https://www.nigms.nih.gov), the [National Institute of General Medical Sciences of the National Institutes of Health](https://nih.gov), and the [National Science Foundation](https://nsf.gov).

<div class="logos">
    <div class="logos-row">
        <a href="https://nih.gov" rel="noopener" target="_blank" title="NIH">
            <img class="zoom" src="/assets/images/about/funding/nih.svg" />
        </a>
        <a href="https://nibib.nih.gov" rel="noopener" target="_blank" title="NIBIB">
            <img class="zoom" src="/assets/images/about/funding/nibib.svg" />
        </a>
        <a href="https://nigms.nih.gov" rel="noopener" target="_blank" title="NIGMS">
            <img class="zoom" src="/assets/images/about/funding/nigms.svg" />
        </a>
        <a href="https://nsf.gov" rel="noopener" target="_blank" title="NSF">
            <img class="zoom" src="/assets/images/about/funding/nsf.svg" />
        </a>
    </div>
</div>
# Citing BioSimulations and BioSimulators

## BioSimulations
Please check back soon for citation information!

## runBioSimulations
Please use the following article to cite runBioSimulations:

Bilal Shaikh, Gnaneswara Marupilla, Mike Wilson, Michael L Blinov, Ion I Moraru & Jonathan R Karr. RunBioSimulations: an extensible web application that simulates a wide range of computational modeling frameworks, algorithms, and formats. Nucleic Acids Research 49(W1): W597-W602 (2021). DOI: [10.1093/nar/gkab411](https://doi.org/10.1093/nar/gkab411)

## BioSimulators
Please check back soon for citation information!
# Format for reports of the results of simulation experiments described with SED-ML

## Overview

The BioSimulators/BioSimulations format for simulation results outlines how SED-ML reports (`sedml:report`) and plots (`sedml:plot2D`, `sedml:plot3D`) should be encoded into [Hierarchical Data Format (HDF) 5](https://www.hdfgroup.org/solutions/hdf5/). These conventions are capable of capturing reports and plots with multiple dimensions, and with data sets that have different shapes and data types, and repeated labels.


## Specifications

Data for reports and plots of simulation results should be saved in HDF5 according to the following conventions:

- Paths of reports and plots: Within the HDF5 file, each report/plot should be saved to a path equal to the combination of (a) the relative location of the parent SED-ML document within the parent COMBINE/OMEX archive and (b) the `id` of the report/plot. For example, a report with id `time_course_results` in a SED-ML file located at `./path/to/experiment.sedml` should be saved to the path `path/to/experiment.sedml/time_course_results` in the HDF5 file.

- Data set shapes: For SED-ML reports, the rows of each HDF5 dataset should correspond to the SED-ML data sets (`sedml:dataSet`) specified in the SED-ML definition of the report (e.g., time symbol, specific model variables). For SED-ML plots, the rows of each HDF5 dataset should correspond to the SED-ML data generators (`sedml:dataGenerator`) specified in the SED-ML definition of the plot (e.g., time symbol, specific model variables).

- Elemental tasks (`sedml:task`):
    - Steady-state simulations (`sedml:steadyState`): The rows of HDF5 data sets should be scalars.
    - One step simulations (`sedml:oneStep`): The rows of HDF5 data sets should be tuples of the start and end points of the simulation.
    - Time course simulations (`sedml:uniformTimeCourse`): The rows of HDF5 data sets should be a vector with length equal to the number of steps of the time course plus one.
    - Simulations of spatial models: The rows of HDF5 data sets should be matrices whose dimensions represent space and time.
    
- Repeated tasks (`sedml:repeatedTask`): The first dimension of each row should represent the iterations of the tasks that produced its values. The second dimension of each data set should represent the individual sub-tasks (`sedml:subTask`) of the task. The results of sub-tasks should be ordered in the same order -- the order of their attributes -- that the sub-tasks were executed. If repeated tasks are nested within repeated tasks, the next dimensions should alternate between representing the iterations and sub-tasks of the nested repeated tasks. The final dimensions of each row should be encoded as above for `sedml:task`. For example, non-spatial time course simulations should have a single additional dimension of length equal to the number of steps of the time course plus one.

!!!note
    If the rows of an HDF5 data set have different shapes, the data sets should be reshaped into a consistent shape by right-padding their values with NaN.

- Metadata for reports: The following metadata should be encoded into attributes of the corresponding HDF5 dataset.

    - Type of the output: The type of the output (`Report`, `Plot2D`, `Plot3D`) should be encoded into the key `_type`.
    - Complete id of the output: The complete id of the output (combination of the location of the parent SED-ML file of the output (`omex-manifest:content/@location`) within its parent COMBINE/OMEX archive and the SED-ML id of the output (`sedml:output/@sedml:id`)) should be encoded into the key `uri`.
    - Id of the output: The SED-ML id of the output (`sedml:output/@sedml:id`) should be encoded into the key `sedmlId`.
    - Name of the output: The name of the output (`sedml:output/@sedml:name`) should be encoded into the key `sedmlName`.
    - Ids of rows (SED-ML data sets or data generators): For reports, the ids of the data sets should be encoded into the key `sedmlDataSetIds`. The value of this key should be an array of the ids of the data sets, in the order in which the data sets were defined in their parent SED-ML document. For plots, the ids of the data generators should be encoded into the key `sedmlDataSetIds`. The value of this key should be an array of the ids of the data generators, in the order in which the data generators were defined in their parent SED-ML document.
    - Names of row (SED-ML data sets or data generators): For reports, the names of the data sets should be encoded into the key `sedmlDataSetNames`. For plots, the names of the data generators should be encoded into the key `sedmlDataSetNames`. The value of this key should be an array of the ids of the data sets, in the order in which the data sets were defined in their parent SED-ML document.
    - Labels of rows (SED-ML data sets or data generators): For reports, the labels of the data sets should be encoded into the key `sedmlDataSetLabels`. For plots, the id of the data generators should be encoded into the key `sedmlDataSetLabels`. The value of this key should be an array of the labels of the data sets, in the order in which the data sets were defined in their parent SED-ML document.
    - Data types of SED-ML data sets/generators: The data types of the data sets (reports) or data generators (plots) should be encoded into the key `sedmlDataSetDataTypes`. The value of this key should be an array of the data types of the data sets/generators, in the order in which the data sets/generators were defined in their parent SED-ML document. The data type of each data set should either be described using a NumPy `dtype` (e.g., `int64`) to indicate a data set whose value is non-null or `__None__` to indicate a data set whose value is `null`.
    - Shapes of SED-ML data sets/generators: The shapes of the data sets (reports) or data generators (plots) should be encoded into the key `sedmlDataSetShapes`. The value of this key should be an array of comma-separated lists of the shapes of the data sets/generators. The shapes of the data sets/generators should be listed in the order in which the data sets/generators were defined in their parent SED-ML document.

- Metadata for SED-ML files: The following metadata should be encoded into attributes of the parent groups of HDF5 datasets which represent SED-ML files and their parent directories within their parent COMBINE archives.

    - Complete id of the COMBINE archive location: The location of each SED-ML file and the location of each parent directory of each SED-ML file with their parent COMBINE archive (`omex-manifest:content/@location`) should be encoded into the keys `uri` and `combineArchiveLocation`.

## Example HDF5 report files

Several example reports are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples).

Below is a graphical illustration of the organization of an HDF5 file for a SED-ML report with id `report-1` defined in a SED-ML file located at `experiment-1/batch-1/simulation-1.sedml` within a COMBINE/OMEX archive.

#### Path of the HDF5 dataset for the SED-ML report
`experiment-1/batch-1/simulation-1.sedml/report-1`

#### HDF5 dataset for the SED-ML report

![example-report](./images/hdf.png)

#### Attributes of the HDF5 dataset for the SED-ML report
- `_type`: `Report`
- `uri`: `experiment-1/batch-1/simulation-1.sedml/report-1`
- `sedmlId`: `report-1`
- `sedmlName`: `Report 1`
- `sedmlDataSetIds`: `time`, `metabolite_a`, `metabolite_b`, `sum_metabolite_a_b`, `ratio_flux_c_d`
- `sedmlDataSetLabels`: `Time`, `Metabolite A`, `Metabolite B`, `Sum of metabolites A and B`, `Flux ratio of reactions C and D`
- `sedmlDataSetDataTypes`: `float64`, `int64`, `int64`, `int64`, `float64`
- `sedmlDataSetShapes`: `14`, `9`, `11`, `9`, `14`

#### Attributes of the HDF5 groups for the SED-ML file and its parent subdirectories
- `experiment-1` HDF5 group for the grandparent directory of the SED-ML file
  - `uri`: `experiment-1`
  - `combineArchiveLocation`: `experiment-1`
- `experiment-1/batch-1` HDF5 group for the parent directory of the SED-ML file
  - `uri`: `experiment-1/batch-1`
  - `combineArchiveLocation`: `experiment-1/batch-1`
- `experiment-1/batch-1/simulation-1.sedml` HDF5 group for the SED-ML file
  - `uri`: `experiment-1/batch-1/simulation-1.sedml`
  - `combineArchiveLocation`: `experiment-1/batch-1/simulation-1.sedml`

## Recommended resources for building reports of simulation results

Below are helpful tools for building reports of simulation results:

- [BioSimulators utils](https://docs.biosimulators.org/Biosimulators_utils/) is a Python library which provides functions for generating reports to the above specifications.
- [h5py](https://www.h5py.org/)  is a high-level Python library for reading and writing HDF5 files.
- [HDF5 libraries](https://www.hdfgroup.org/downloads/hdf5) for C, C++, and Java.
# Conventions

BioSimulators and BioSimulations follow a set of conventions designed to foster consistent representation and interpretation of simulation projects across modeling frameworks, simulation algorithms, modeling formats, and simulation tools. In particular, these conventions ensure that simulation tools produce consistent outputs (reports and plots at consistent paths in consistent formats) from consistent inputs (descriptions of models and simulations). These conventions also help investigators (a) communicate the semantic meaning of simulation projects and the capabilities of simulation tools, (b) find simulation projects that represent specific biology and simulation tools that have specific capabilities, and (c) reproduce and reuse simulation projects with alternative software tools.

These conventions encompass multiple formats and guidelines:

- [Format for describing the specifications of simulation tools](./simulator-capabilities.md): This format enables investigators to describe the modeling frameworks, simulation algorithms, and modeling formats supported by a simulation tool. The format also enables investigators to describe the parameters of each algorithm, their data types, and their default values.

- [Conventions for command-line applications and APIs for simulation tools](./simulator-interfaces.md): These conventions ensure that simulation tools (a) support consistent syntax for executing simulations and (b) produce outputs at consistent locations in consistent formats.

- [Conventions for Docker images for simulation tools](./simulator-images.md): These conventions ensure that the entry points of containerized simulation tools support consistent syntax for executing simulations and that containerized simulation tools provide consistent metadata.

- [Conventions for simulation experiments with SED-ML](./simulation-experiments.md): These conventions ensure that the community consistently encodes simulation experiments into SED-ML. This includes conventions for targets for implicit elements of models, or symbols, which are not directly defined in models (e.g., reduced costs of FBA reactions, shadow prices of FBA species). This also delineates how to encode the values of model attribute changes and algorithm parameters into SED-ML, including encoding enumerated values, lists, dictionaries, and other data structures.

- [Format for reports of simulation results](./simulation-run-reports.md) (`sedml:report`): This format ensures that simulation tools produce reports and plots in a consistent format (e.g., [HDF5](https://www.hdfgroup.org/solutions/hdf5/)) with consistent shapes (e.g., rows: data set (`sedml:dataSet`), columns: simulation step) that can be consistently visualized and interpreted.

- [Guidelines for data visualizations with Vega](./simulation-run-visualizations.md): [Vega](https://vega.github.io/vega/) is a powerful format for describing interactive, two-dimensional, data visualizations. Vega makes visualizations re-usable by separately capturing visual marks and how they should be painted with data. These guidelines outline how to use Vega to visualize the results of simulation experiments captured by SED-ML reports.

- [Guidelines for using the OMEX Metadata format to annotate the meaning and provenance of simulation projects](./simulation-project-metadata.md) (COMBINE/OMEX archives): We recommend using the OMEX Metadata RDF-XML format to annotate the meaning, provenance, and credibility of simulation projects. These guidelines recommend specific predicates and objects for annotating simulation projects and their components.

- [Format for logs of the execution of simulation projects](./simulation-run-logs.md) (COMBINE/OMEX archives): This format enables simulation tools to communicate the status and outcome of the execution of simulation projects. This format enables simulation tools to communicate information such as the following:

    - The status and outcome of the project and each individual SED-ML document, task, and output (e.g., succeeded, failed, skipped, queued).
    - The function and arguments used to execute each individual simulation step.
    - The standard output/error produced by the execution of the project and each SED-ML document, task, and output.
    - The reason why each SED-ML document, task or output could not be executed.
     
# Format for logs of the execution of simulation experiments captured in COMBINE/OMEX archives

## Overview

This format enables simulation tools to communicate information about the status and outcome of the execution of SED-ML files in COMBINE/OMEX archives. The format can capture the following information:

- The status and outcome of the COMBINE archive and each SED-ML document, task, report, plot, data set, curve and surface (e.g., queued, running, succeeded, failed).
- Information about simulation functions that were executed and the arguments that were used.
- The standard output/error produced from executing the COMBINE archive and each SED-ML document, task, and output.
- The duration of the execution of the COMBINE archive and each SED-ML document, task, and output.
- The reason for each skipped or failed SED-ML document, task or output.

## Specifications

The schema for the format is available in [JSON schema](https://api.biosimulations.org/schema/CombineArchiveLog.json) and [Open API](https://api.biosimulations.org/openapi.json) formats. Documentation for the schema is available at [https://api.biosimulations.org/](https://api.biosimulations.org/).

Simulators are encouraged to log the execution of each individual SED-ML element, including each task, report, plot, data set, curve and surface.

At the same time, the format provides simulation tools flexibility to log their execution at whatever level of granularity is possible. Below are several possible levels of granularity.

- Individual task and output elements (e.g., data sets, curves, surfaces)
- Individual tasks
- Individual SED-ML documents
- Entire COMBINE/OMEX archives

!!!note 
    The output attribute for the log of each COMBINE/OMEX archive, SED-ML document, task, and output can include [ANSI escape codes](https://en.wikipedia.org/wiki/ANSI_escape_code) for color and other terminal formatting.


## Example logs

Below is an example of an element-level log of a COMBINE/OMEX archive that involves a single SED-ML document with two tasks at the beginning of its execution.

```json
{
  "status": "QUEUED",
  "exception": null,
  "skipReason": null,
  "output": null,
  "duration": null,
  "sedDocuments": [
    {
      "location": "doc_1.sedml",
      "status": "QUEUED",
      "exception": null,
      "skipReason": null,
      "output": null,
      "duration": null,
      "tasks": [
        {
          "id": "task_1_ss",
          "status": "QUEUED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": null,
          "algorithm": null,
          "simulatorDetails": null
        },
        {
          "id": "task_2_time_course",
          "status": "QUEUED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": null,
          "algorithm": null,
          "simulatorDetails": null
        }
      ],
      "outputs": [
        {
          "id": "report_1",
          "status": "QUEUED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": null,
          "dataSets": [
            {
              "id": "dataset_1",
              "status": "QUEUED"
            },
            {
              "id": "dataset_2",
              "status": "QUEUED"
            }
          ]
        },
        {
          "id": "plot_1",
          "status": "QUEUED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": null,
          "curves": [
            {
              "id": "curve_1",
              "status": "QUEUED"
            }
          ]
        }
      ]
    }
  ]
}
```

Below is an example of an element-level log of the final state of the successful execution of the same COMBINE/OMEX archive.

```json

{
  "status": "SUCCEEDED",
  "exception": null,
  "skipReason": null,
  "output": null,
  "duration": 6,
  "sedDocuments": [
    {
      "location": "doc_1.sedml",
      "status": "SUCCEEDED",
      "exception": null,
      "skipReason": null,
      "output": null,
      "duration": 5,
      "tasks": [
        {
          "id": "task_1_ss",
          "status": "SUCCEEDED",
          "exception": null,
          "skipReason": null,
          "output": "Reading model ... done\nInitializing simulation ... done\nExecuting simulation ... done\n",
          "duration": 2,
          "algorithm": null,
          "simulatorDetails": null
        },
        {
          "id": "task_2_time_course",
          "status": "SUCCEEDED",
          "exception": null,
          "skipReason": null,
          "output": "Reading model ... done\nInitializing simulation ... done\nExecuting simulation ... done\n",
          "duration": 1,
          "algorithm": null,
          "simulatorDetails": null
        }
      ],
      "outputs": [
        {
          "id": "report_1",
          "status": "SUCCEEDED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": 0.1,
          "dataSets": [
            {
              "id": "dataset_1",
              "status": "SUCCEEDED"
            },
            {
              "id": "dataset_2",
              "status": "SUCCEEDED"
            }
          ]
        },
        {
          "id": "plot_1",
          "status": "SUCCEEDED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": 0.01,
          "curves": [
            {
              "id": "curve_1",
              "status": "SUCCEEDED"
            }
          ]
        }
      ]
    }
  ]
}
```

Below is an example of an element-level log of the final state of the failed execution of the same COMBINE/OMEX archive.

```json
{
  "status": "FAILED",
  "exception": null,
  "skipReason": null,
  "output": null,
  "duration": 6,
  "sedDocuments": [
    {
      "location": "doc_1.sedml",
      "status": "FAILED",
      "exception": null,
      "skipReason": null,
      "output": null,
      "duration": 5,
      "tasks": [
        {
          "id": "task_1_ss",
          "status": "SUCCEEDED",
          "exception": null,
          "skipReason": null,
          "output": "Reading model ... done\nInitializing simulation ... done\nExecuting simulation ... done\n",
          "duration": 2,
          "algorithm": null,
          "simulatorDetails": null
        },
        {
          "id": "task_2_time_course",
          "status": "FAILED",
          "exception": {
            "type": "FileNotFoundError",
            "message": "Model `model2.xml` does not exist."
          },
          "skipReason": null,
          "output": null,
          "duration": 1,
          "algorithm": null,
          "simulatorDetails": null
        }
      ],
      "outputs": [
        {
          "id": "report_1",
          "status": "SUCCEEDED",
          "exception": null,
          "skipReason": null,
          "output": null,
          "duration": 0.1,
          "dataSets": [
            {
              "id": "dataset_1",
              "status": "SUCCEEDED"
            },
            {
              "id": "dataset_2",
              "status": "SUCCEEDED"
            }
          ]
        },
        {
          "id": "plot_1",
          "status": "SKIPPED",
          "exception": null,
          "skipReason": {
            "type": "2DPlotNotImplemented",
            "message": "Output skipped because the simulator cannot generate plots."
          },
          "output": null,
          "duration": 0.01,
          "curves": [
            {
              "id": "curve_1",
              "status": "SKIPPED"
            }
          ]
        }
      ]
    }
  ]
}
```

Below is an example of a SED-ML document-level log of the final state of the failed execution of the same COMBINE/OMEX archive.

```json

{
  "status": "FAILED",
  "exception": null,
  "skipReason": null,
  "output": null,
  "duration": 6,
  "sedDocuments": [
    {
      "location": "doc_1.sedml",
      "status": "FAILED",
      "exception": {
        "type": "FileNotFoundError",
        "message": "Model `model2.xml` does not exist."
      },
      "skipReason": null,
      "output": "Reading model ... done\nInitializing simulation ... done\nExecuting simulation ... done\n",
      "duration": 5,
      "tasks": null,
      "outputs": null
    }
  ]
}

```

## Events which should trigger simulators to update the status of COMBINE/OMEX archives

Simulation tools are encouraged to flush logs of the execution of COMBINE/OMEX archives upon each of the following events:

- Start of the execution of the COMBINE/OMEX archive. After this event, the status of the archive should be `RUNNING`, the status of each SED-ML document and SED-ML element that the simulation tool is capable of executing should be `QUEUED`, and the status of each SED-ML document and SED-ML element that the simulation tool is not capable fo executing should be `SKIPPED`.
- Start of the execution of each SED-ML document. This event should change the status of the document to `RUNNING`.
- Start of the execution of each SED-ML task. This event should change the status of the task to `RUNNING`.
- End of the execution of each SED-ML task. If the task succeeded, its status should be changed to `SUCCEEDED`. In addition, all data sets, curves, and surfaces which can now be generated should be generated, and their status should be changed to `SUCCEEDED`. The status of the parent reports and plots should be set to `SUCCEEDED` if all of their data sets, curves, and surfaces have been generated, or `RUNNING` if some of their data sets, curves, and surfaces cannot yet be generated because they depend on tasks which have not yet been executed. If the task failed, its status should be changed to `FAILED`, and the status of all data sets, curves, surfaces, reports, and plots which depend on the task should be changed to `FAILED`.
- End of the execution of each SED-ML document. If all of the task and outputs in the document succeeded, the document's status should be changed to `SUCCEEDED`. Otherwise, the document's status should be changed to `FAILED`.

By the end of the execution of a COMBINE/OMEX archive, the status of each SED-ML document, task, report, plot, data set, curve, and surface should be one of `SUCCEEDED`, `SKIPPED`, or `FAILED`.


## Recommended resources for building logs of the execution of COMBINE/OMEX archives

Below are helpful tools for building logs of the execution of COMBINE/OMEX archives:

- [BioSimulators utils](https://docs.biosimulators.org/Biosimulators_utils/) is a Python library which provides functions for generating reports to the above specifications.
- [capturer](https://pypi.org/project/capturer/) is a Python library for capturing standard output and standard error streams.
# Standard command-line applications and Python APIs for biosimulation tools

## Overview

The BioSimulators conventions for command-line applications and Python APIs for biosimulation tools are sets of requirements for the syntax and semantics of the inputs and outputs of biosimulation software tools. The conventions ensure that simulation tools can be consistently executed with the same input arguments (a path to a COMBINE/OMEX archive that defines models and simulations, and a path to save the outputs of the simulation experiments defined in the archive), and that the simulation tools produce consistent outputs (reports and plots at consistent paths in consistent formats).

We encourage developers to provide two interfaces for two purposes:

- Command-line entrypoints to Docker images: These interfaces provide developers a wide degree of flexibility (e.g., to use different programming languages and dependencies), make it easy for BioSimulators to archive each version of each tool, and make it easy for investigators to use tools to execute COMBINE/OMEX archives. However, these entrypoints are not optimally efficient for all use cases, and they provide investigators limited flexibility.
- Python APIs: APIs provide investigators additional flexibility, such as to develop higher-performance simulation services or to co-simulate models using multiple tools. However, for optimal efficiency, these APIs require developers to use a specific programming language and to package their tools for easy installation or provide straightforward installation instructions. Furthermore, these APIs are typically more complex for investigators to install and use.

In most cases, we recommend that developers create BioSimulators-compliant interfaces in three steps:

1. Use [BioSimulators-utils](https://github.com/biosimulators/Biosimulators_utils) to develop a compliant Python API.
2. Use `biosimulators_utils.simulator.cli.build_cli` to construct a command-line application from their Python API.
3. Construct a Docker image by attaching their command-line application to the entrypoint of an image.
4. For validation, simulation tools MUST provide Docker images whose entrypoints are BioSimulators-compliant command-line applications. Simulation tools are also OPTIONALLY encouraged to provide BioSimulators-compliant APIs.

## Convention for command-line applications

Simulation tools should support the following command-line arguments:

- `-i, --archive`: A path to a COMBINE/OMEX archive which contains descriptions of one or more simulation tasks.
- `-o, --out-dir`: The path where the outputs (reports and plots) of the simulation tasks should be saved.
- `-h, --help`: An optional argument that instructs the command-line application to print help information about itself.
- `-v, --version`: An optional argument that instructs the command-line application to report its version and the versions of any critical dependencies.

As an example, below is the documentation for the command-line application for the [tellurium](https://biosimulators.org/simulators/tellurium)  biochemical simulation program.

```bash
usage: tellurium [-h] [-d] [-q] -i ARCHIVE [-o OUT_DIR] [-v]

BioSimulators-compliant command-line interface to the tellurium simulation program <http://tellurium.analogmachine.org>.

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           full application debug mode
  -q, --quiet           suppress all console output
  -i ARCHIVE, --archive ARCHIVE
                        Path to a COMBINE/OMEX archive file which contains one or more SED-ML-
                        encoded simulation experiments
  -o OUT_DIR, --out-dir OUT_DIR
                        Directory to save outputs
  -v, --version         show program's version number and exit
```

## Conventions for Python APIs

Python APIs should be Python modules which provide the following attributes and methods:

- `__version__ (str)`: Version of the API (e.g., `0.1.0`).
- `get_simulator_version (() -> str)`: A function which get the versions of the underlying simulation tool (e.g., `1.0.0`).
- `preprocess_sed_task ((task: Task, variables: list[Variable]) -> Any)`: A function which preprocesses the information required to execute a SED-ML task. This method enables users to efficiently execute multiple simulation steps (using `exec_sed_task`) without unnecessary duplicate computations, such as to import models from files, validate models, and identify suitable algorithms.

- `exec_sed_task ((task: Task, variables: list[Variable], preprocessed_task:Any=None) -> Tuple[VariableResults, TaskLog])`: A function which executes a single SED-ML task involving a single simulation of a single model and returns the predicted value of each requested output variable and a log of the execution of the simulation.

- `exec_sed_doc ((doc: SedDocument, variables: list[Variable]) -> Tuple[ReportResults, SedDocumentLog])`: A function which executes a single SED-ML document involving one or more simulations of one or more models, returns data for each SED-ML report and plot, and returns a log of the execution of the simulations, reports, and plots.

- `exec_sedml_docs_in_combine_archive ((archive_filename: str, out_dir: str, return_results: bool = False) -> Tuple[SedDocumentResults, CombineArchiveLog])`: A function which executes all of the tasks in all of the SED-ML files in a COMBINE/OMEX archive, exports all of the requested reports and plots, optionally returns the result of each output, and returns a log of the execution of the archive.

More information about the required method signatures is available in the [simulator template](https://github.com/biosimulators/Biosimulators_simulator_template/blob/dev/my_simulator/__init__.py).

## Execution behavior

Both command-line interfaces and Python APIs should also support the following conventions:

- Support at least one "standard" modeling language. Interfaces to simulation tools should support at least a subset of at least one community modeling language such as BNGL, CellML, Kappa, NeuroML/LEMS, HOC, pharmML, SBML, Smoldyn, VCML, or XPP.
- Support the SED-ML format for describing simulation experiments. Interfaces to simulation tools should support at least the following features of SED-ML:
    - Models and model attribute changes: `sedml:model`, `sedml:changeAttribute`
    - At least one of steady-state, one step, or timecourse simulations: `sedml:steadyState`, `sedml:oneStep`, or `sedml:uniformTimeCourse`.
    - Tasks to execute a single simulation of a single model: `sedml:task`.
    - Algorithms and their parameters: `sedml:algorithm`, `sedml:algorithmParameter`.
    - Data generators for individual variables: `sedml:dataGenerator`
    - Report outputs: `sedml:report`.
- Use KiSAO to describe simulation algorithms and their parameters. Interfaces to simulation tools should use KiSAO terms to indicate specific algorithms and algorithm parameters.
- Support the COMBINE/OMEX format for collections of models and simulations. Interfaces to simulation tools should support the full COMBINE/OMEX specification.
- Execute all tasks described in the master file of the COMBINE/OMEX archive: When a COMBINE/OMEX archive file contains a master SED-ML document, simulation tools should execute all tasks defined inside the master file. When an archive doesn't contain a master SED-ML file, simulation tools should execute all of the tasks defined in each SED-ML file in the archive.
- Support the BioSimulators format for reports of simulation results: Interfaces to simulation tools should save reports in the [BioSimulators HDF5 format](./simulation-run-reports.md) for simulation data from reports and plots. Within the HDF5 file, each report and plot should be saved to a path equal to the combination of the relative path of its parent SED-ML file within the COMBINE/OMEX archive and the id of the report.
- Save plots in Portable Document Format (PDF) bundled into a zip archive. Within the zip archive, each plot should be saved to a path equal to the combination of the relative path of its parent SED-ML file within the COMBINE/OMEX archive, the id of the plot, and the extension `.pdf`.
- Save simulation outputs to standard file paths: Data for reports and plots should be saved to `{ out-dir }/reports.h5`. Plots should be saved to `{ out-dir }/plots.zip`.

## Environment variables
To further support consistent execution of simulations with other simulation tools, command-line interfaces and Python APIs are also encouraged to implement the following environment variables. The Dockerfiles for simulation tools should use the `ENV` directive to indicate the variables they support and their default values.

- `ALGORITHM_SUBSTITUTION_POLICY`: This environment variable enables users to control if and how the simulator substitutes algorithms with other mathematically-equivalent or similar algorithms.

    We recognize the increasing levels of substitution listed below. Simulation tools are encouraged to use `SIMILAR_VARIABLES` as the default value for `ALGORITHM_SUBSTITUTION_POLICY`.

    A recommended matrix of algorithm substitutions is available from the [KiSAO](https://github.com/SED-ML/KiSAO/blob/dev/libkisao/python/docs/algorithm-substitutability.csv) ontology. A [Python package](https://github.com/SED-ML/KiSAO/blob/dev/libkisao/python/) is also available for implementing this matrix.

    When alternate algorithms are substituted, we recommend that simulation tools ignore SED-ML algorithm parameters as algorithm parameters can have different meanings in the context of different algorithms.

    For algorithm substitution level `NONE`, we recommend that simulation tools raise errors for unsupported algorithm parameters and unsupported values of algorithm parameters. For higher substitution levels, we recommend that simulation tools skip unsupported parameters and unsupported values and raise warnings when parameters are skipped.

    0. `NONE`: Tools should strictly interpret SED-ML simulations, and raise errors on the execution of SED-ML files that involve unsupported algorithms. For example, a simulation tool that only supports the Stochastic Simulation Algorithm (SSA, KISAO_0000029) should raise errors on SED-ML files that request simulations with the Next Reaction Method (NRM, KISAO_0000027). In many cases, this level will effectively constrain the execution of a SED-ML document to a specific implementation of an algorithm by a specific simulation tool.

    1. `SAME_METHOD`: Algorithms can be substituted with different realizations of the same method. For example, GLPK's implementation of the Simplex method could be substituted with SciPy's implementation.

    2. `SAME_MATH`: Tools should execute simulations with alternative mathematically-equivalent algorithms, and raise errors on the execution of SED-ML files which request algorithms that are mathematically distinct from those implemented by the tool. When tools execute alternative mathematically-equivalent algorithms, they should issue warnings to this effect. For example, a simulation tool that only supports SSA should execute simulations that request NRM with a warning, and raise an error on SED-ML files that request the tau-leaping method (KISAO_0000039).

    3. `SIMILAR_APPROXIMATIONS`: Algorithms can be substituted with others that make similar approximations to the same mathematics. For example, CVODE could be substituted with LSODA or the Fehlberg method. Tau leaping could be substituted with partitioned tau leaping.

    4. `DISTINCT_APPROXIMATIONS`: Algorithms can be substituted with others that make distinct approximations to the same math. For example, SSA could be substituted with tau leaping or the Pahle hybrid method.

    5. `DISTINCT_SCALES`: Algorithms can be substituted with others that make distinct approximations to the same math that substantially differ in their scales. For example, SSA could be substituted with CVODE.

    6. `SAME_VARIABLES`: Algorithms that predict the same output variables can be substituted. For example, FBA could be substituted with parsimonious FBA.

    7. `SIMILAR_VARIABLES` (recommended default): Algorithms that predict similar output variables can be substituted. For example, FBA could be substituted with geometric FBA.

    8. `SAME_FRAMEWORK`: Tools should execute simulations with alternative algorithms, including algorithms that are not mathematically equivalent, and issue warnings when alternative algorithms are executed.

    9. `ANY`: Tools can execute simulations with any alternative algorithm. Note, switching to any other algorithm can substantially change the interpretation of a simulation (e.g., switching SSA to CVODE loses all information about the variance of a simulation).

- `VERBOSE`: Indicates whether a simulator should display detailed information about the execution of each task.

    We recognize the following values.

    - `0`, `false` (any case)
    - `1`, `true` (any case)

## Execution of modeling projects encoded as COMBINE/OMEX archives

To ensure consistent execution of simulation experiments, command-line applications and Python APIs should adopt the conventions described below for the execution of COMBINE/OMEX archives.

- Identification of SED-ML files. SED-ML files should be identified as `omex:content` whose `format` attribute starts with `http://identifiers.org/combine.specifications/sed-ml`.
 
- Preferential execution of "master" files. The OMEX format supports the declaration of a single "master" file (`omex:content[@master='true']`). When a COMBINE/OMEX archive contains a single master file, simulation tools should only execute this file. Note, if the master file is not a SED-ML document, then no simulations should be executed. When a COMBINE/OMEX archive doesn't have a master file, all SED-ML documents should be executed.

## Execution of simulation experiments encoded in SED-ML

To ensure consistent execution of simulation experiments, command-line applications and Python APIs should adopt the conventions described below for the execution of SED-ML files.

- Substitution of alternative simulation algorithms.

    Because no simulation tool implements every simulation algorithm, simulation tools are encouraged to execute SED-ML simulations with alternative algorithms (close KiSAO terms) when the tool does not support the requested algorithm (`sedml:algorithm/@sedml:kisaoID`). For example, a tool which only implements the Stochastic Simulation Algorithm (SSA, `KISAO_0000029`) could choose to execute simulations that request the Next Reaction Method (NRM, `KISAO_0000027`), a mathematically-equivalent method, with SSA.

    Simulation tools are encouraged to use the KiSAO ontology to systematically identify related simulation algorithms.

    When a tool uses an alternative algorithm, the tool should issue a warning message to the user that indicates that an alternative algorithm was used.

    Tools which choose to execute alternative algorithms should support the `ALGORITHM_SUBSTITUTION_POLICY` environment variable (see [above](#environment-variables)).

## Recommended resources for implementing command-line applications and APIs

Below are helpful tools for implementing command-line applications and Python APIs to the above specifications:

- [BioSimulators utils](https://docs.biosimulators.org/Biosimulators_utils/) is a Python library which provides functions implementing command-line applications to the above specifications, as well as functions for interpreting COMBINE/OMEX archives and SED-ML files, generating tables and plots of simulation plots, and logging the execution of COMBINE/OMEX archives. BioSimulators utils provides high-level access to some of the lower-level libraries listed below.
- [libCOMBINE](https://github.com/sbmlteam/libCombine) is a library for creating and unpacking COMBINE/OMEX archives. libCOMBINE provides bindings for several languages.
- [libSED-ML](https://github.com/fbergmann/libSEDML) is a library for serializing and deserializing SED-ML documents to and from XML files. libSED-ML provides bindings for several languages.
- [libOmexMeta](https://github.com/sys-bio/libOmexMeta) is a library for reading and querying OMEX Metadata files. libOmexMeta provides bindings for several languages.
- [argparse](https://docs.python.org/3/library/argparse.html) is a Python module for implementing command-line applications.
- [Cement](https://builtoncement.com/) is a higher-level Python library for implementing more complex command-line applications.
# Format for the specification of BioSimulation tools

## Overview
The BioSimulators format for the specifications of a simulation tool is a JSON schema for describing the modeling frameworks (e.g., logical, constraint-based), simulation algorithms (e.g., CVODE, SSA), and modeling formats (e.g., CellML, COMBINE/OMEX, SBML, SED-ML) that a simulation tool supports, as well as the parameters of each algorithm (e.g., random number generator seed), their data types, and their allowed and default values.

The format can also capture a metadata about each simulation tool including its

- Name;
- Version;
- Description;
- URL for a standardized Docker image for the tool;
- URL for documentation about the tool;
- Citations for the tool;
- Citations for each algorithm supported by the tool;
- Other identifiers for the tool such as for bio.tools, CRAN, or PyPI;
- License for the tool;
- The authors of the tool; and
- Dates when the tool was submitted to the BioSimulators registry and when it was last updated.

## Schema

The schema for the format is available in [JSON Schema](https://api.biosimulators.org/schema/Simulator.json) and [Open API](https://api.biosimulators.org/openapi.json) formats. Documentation for the schema is available at [https://api.biosimulators.org/](https://api.biosimulators.org/).

The schema utilizes several ontologies:

- Funding agencies: [Funder Registry](https://www.crossref.org/services/funder-registry/) terms such as the National Science Foundation ([10.13039/100000001](http://doi.org/10.13039/100000001)) and the National Institutes of Health ([10.13039/100000002](http://doi.org/10.13039/100000002)).
 
- Licenses: [SPDX](https://spdx.org/) terms such as GNU General Public License v3.0 or later ([GPL-3.0-or-later](https://spdx.org/licenses/GPL-3.0-or-later)).

- Modeling formats: [EDAM](https://edamontology.org/) terms such as BNGL ([format_3972](https://www.ebi.ac.uk/ols/ontologies/edam/terms?iri=http%3A%2F%2Fedamontology.org%2Fformat_3972)), CellML ([format_3240](https://www.ebi.ac.uk/ols/ontologies/edam/terms?iri=http%3A%2F%2Fedamontology.org%2Fformat_3240)), SBML ([format_2585](https://www.ebi.ac.uk/ols/ontologies/edam/terms?iri=http%3A%2F%2Fedamontology.org%2Fformat_2585)), and SED-ML ([format_3685](https://www.ebi.ac.uk/ols/ontologies/edam/terms?iri=http%3A%2F%2Fedamontology.org%2Fformat_3685)).

- Modeling frameworks: [SBO](https://www.ebi.ac.uk/sbo/) terms such as flux balance analysis framework ([SBO:0000624](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000624)) and non-spatial continuous kinetic framework ([SBO:0000293](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000293)).

- Programming languages: [Linguist](https://github.com/github/linguist) terms such as C++, Java, JavaScript, Python, R.

- Simulation algorithms: [KiSAO](http://co.mbine.org/standards/kisao) terms such as CVODE ([KISAO:0000019](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000019)), FBA ([KISAO:0000437](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000437)), and SSA ([KISAO:0000029](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000029)).

## Example

As an example, below is the specification for [tellurium](https://biosimulators.org/simulators/tellurium/):

```json
{
  "id": "tellurium",
  "name": "tellurium",
  "version": "2.2.1",
  "description": "tellurium is a Python-based environment for model building, simulation, and analysis that facilitates reproducibility of models in systems and synthetic biology.",
  "urls": [
    {
      "type": "Home page",
      "url": "http://tellurium.analogmachine.org/"
    },
    {
      "type": "Discussion forum",
      "url": "https://groups.google.com/d/forum/tellurium-discuss"
    },
    {
      "type": "Tutorial",
      "url": "https://tellurium.readthedocs.io/en/latest/quickstart.html"
    },
    {
      "type": "Installation instructions",
      "url": "https://github.com/sys-bio/tellurium#installation-instructions"
    },
    {
      "type": "Documentation",
      "url": "https://tellurium.readthedocs.io/"
    },
    {
      "type": "Source repository",
      "url": "https://github.com/sys-bio/tellurium"
    },
    {
      "type": "Guide to contributing",
      "url": "https://github.com/sys-bio/tellurium/blob/develop/CONTRIBUTING.rst"
    },
    {
      "type": "Issue tracker",
      "url": "https://github.com/sys-bio/tellurium/issues"
    },
    {
      "type": "License",
      "url": "https://github.com/sys-bio/tellurium/blob/develop/LICENSE.txt"
    },
    {
      "type": "Release notes",
      "url": "http://tellurium.analogmachine.org/news/"
    }
  ],
  "image": {
    "url": "ghcr.io/biosimulators/biosimulators_tellurium/tellurium:2.2.1",
    "digest": "sha256:5d1595553608436a2a343f8ab7e650798ef5ba5dab007b9fe31cd342bf18ec81",
    "format": {
      "namespace": "EDAM",
      "id": "format_3973",
      "version": "1.2.0",
      "supportedFeatures": []
    },
    "operatingSystemType": "Linux"
  },
  "cli": {
    "packageRepository": "PyPI",
    "package": "biosimulators-tellurium",
    "command": "biosimulators-tellurium",
    "installationInstructions": "https://docs.biosimulators.org/biosimulators_tellurium/installation.html"
  },
  "pythonApi": {
    "package": "biosimulators-tellurium",
    "module": "biosimulators_tellurium",
    "installationInstructions": "https://docs.biosimulators.org/biosimulators_tellurium/installation.html"
  },
  "authors": [
    {
      "firstName": "Jayit",
      "lastName": "Biswas",
      "identifiers": []
    },
    {
      "firstName": "Kiri",
      "lastName": "Choi",
      "identifiers": [
        {
          "namespace": "orcid",
          "id": "0000-0002-0156-8410",
          "url": "https://orcid.org/0000-0002-0156-8410"
        }
      ]
    },
    {
      "firstName": "Wilbert",
      "lastName": "Copeland",
      "identifiers": []
    },
    {
      "firstName": "Caroline",
      "lastName": "Cannistra",
      "identifiers": []
    },
    {
      "firstName": "Alex",
      "lastName": "Darling",
      "identifiers": []
    },
    {
      "firstName": "Nasir",
      "lastName": "Elmi",
      "identifiers": []
    },
    {
      "firstName": "Michal",
      "lastName": "Galdzicki",
      "identifiers": [
        {
          "namespace": "orcid",
          "id": "0000-0002-8392-8183",
          "url": "https://orcid.org/0000-0002-8392-8183"
        }
      ]
    },
    {
      "firstName": "Stanley",
      "lastName": "Gu",
      "identifiers": []
    },
    {
      "firstName": "Totte",
      "lastName": "Karlsson",
      "identifiers": []
    },
    {
      "firstName": "Matthias",
      "lastName": "König",
      "identifiers": [
        {
          "namespace": "orcid",
          "id": "0000-0003-1725-179X",
          "url": "https://orcid.org/0000-0003-1725-179X"
        }
      ]
    },
    {
      "firstName": "J",
      "middleName": "Kyle",
      "lastName": "Medley",
      "identifiers": [
        {
          "namespace": "orcid",
          "id": "0000-0002-9135-0844",
          "url": "https://orcid.org/0000-0002-9135-0844"
        }
      ]
    },
    {
      "firstName": "Herbert",
      "middleName": "M.",
      "lastName": "Sauro",
      "identifiers": [
        {
          "namespace": "orcid",
          "id": "0000-0002-3659-6817",
          "url": "https://orcid.org/0000-0002-3659-6817"
        }
      ]
    },
    {
      "firstName": "Andy",
      "lastName": "Somogyi",
      "identifiers": []
    },
    {
      "firstName": "Lucian",
      "lastName": "Smith",
      "identifiers": [
        {
          "namespace": "orcid",
          "id": "0000-0001-7002-6386",
          "url": "https://orcid.org/0000-0001-7002-6386"
        }
      ]
    },
    {
      "firstName": "Kaylene",
      "lastName": "Stocking",
      "identifiers": []
    }
  ],
  "references": {
    "identifiers": [
      {
        "namespace": "pypi",
        "id": "tellurium",
        "url": "https://pypi.org/project/tellurium/"
      },
      {
        "namespace": "pypi",
        "id": "biosimulators-tellurium",
        "url": "https://pypi.org/project/biosimulators-tellurium/"
      },
      {
        "namespace": "nanohub.resource",
        "id": "tellurium",
        "url": "https://nanohub.org/resources/tellurium"
      }
    ],
    "citations": [
      {
        "title": "tellurium: an extensible Python-based modeling environment for systems and synthetic biology",
        "authors": "Kiri Choi, J. Kyle Medley, Matthias König, Kaylene Stocking, Lucian Smith, Stanley Gua & Herbert M. Sauro",
        "journal": "BioSystems",
        "volume": "171",
        "pages": "74-79",
        "year": 2018,
        "identifiers": [
          {
            "namespace": "doi",
            "id": "10.1016/j.biosystems.2018.07.006",
            "url": "https://doi.org/10.1016/j.biosystems.2018.07.006"
          }
        ]
      }
    ]
  },
  "license": {
    "namespace": "SPDX",
    "id": "Apache-2.0"
  },
  "algorithms": [
    {
      "id": "cvode",
      "name": "C-language Variable-coefficient Ordinary Differential Equation solver",
      "kisaoId": {
        "namespace": "KISAO",
        "id": "KISAO_0000019"
      },
      "modelingFrameworks": [
        {
          "namespace": "SBO",
          "id": "SBO_0000293"
        }
      ],
      "modelFormats": [
        {
          "namespace": "EDAM",
          "id": "format_2585",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "modelChangePatterns": [
        {
          "name": "Change component attributes",
          "types": [
            "SedAttributeModelChange",
            "SedComputeAttributeChangeModelChange",
            "SedSetValueAttributeModelChange"
          ],
          "target": {
            "value": "//*/@*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Add components",
          "types": [
            "SedAddXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Remove components",
          "types": [
            "SedRemoveXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Change components",
          "types": [
            "SedChangeXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        }
      ],
      "simulationFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3685",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "simulationTypes": [
        "SedUniformTimeCourseSimulation"
      ],
      "archiveFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3686",
          "version": "1",
          "supportedFeatures": []
        }
      ],
      "citations": [
        {
          "title": "CVODE, a stiff/nonstiff ODE solver in C",
          "authors": "Scott D. Cohen, Alan C. Hindmarsh & Paul F. Dubois",
          "journal": "Computers in Physics",
          "volume": "10",
          "issue": "2",
          "pages": "138-143",
          "year": 1996,
          "identifiers": [
            {
              "namespace": "doi",
              "id": "10.1063/1.4822377",
              "url": "https://doi.org/10.1063/1.4822377"
            }
          ]
        }
      ],
      "parameters": [
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000209"
          },
          "id": "relative_tolerance",
          "name": "Relative tolerance",
          "type": "float",
          "value": "0.000001",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000211"
          },
          "id": "absolute_tolerance",
          "name": "Absolute tolerance",
          "type": "float",
          "value": "1e-12",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000220"
          },
          "id": "maximum_bdf_order",
          "name": "Maximum Backward Differentiation Formula (BDF) order",
          "type": "integer",
          "value": "5",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000219"
          },
          "id": "maximum_adams_order",
          "name": "Maximum Adams order",
          "type": "integer",
          "value": "12",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000415"
          },
          "id": "maximum_num_steps",
          "name": "Maximum number of steps",
          "type": "integer",
          "value": "20000",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000467"
          },
          "id": "maximum_time_step",
          "name": "Maximum time step",
          "type": "float",
          "value": null,
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000485"
          },
          "id": "minimum_time_step",
          "name": "Minimum time step",
          "type": "float",
          "value": null,
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000559"
          },
          "id": "initial_time_step",
          "name": "Initial time step",
          "type": "float",
          "value": null,
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000671"
          },
          "id": "stiff",
          "name": "Stiff",
          "type": "boolean",
          "value": "true",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000670"
          },
          "id": "multiple_steps",
          "name": "Multiple steps",
          "type": "boolean",
          "value": "false",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        }
      ],
      "outputDimensions": [
        {
          "namespace": "SIO",
          "id": "SIO_000418"
        }
      ],
      "outputVariablePatterns": [
        {
          "name": "time",
          "symbol": {
            "value": "time",
            "namespace": "urn:sedml:symbol"
          }
        },
        {
          "name": "species concentrations",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species",
            "grammar": "XPath"
          }
        },
        {
          "name": "parameter values",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter",
            "grammar": "XPath"
          }
        },
        {
          "name": "reaction fluxes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction",
            "grammar": "XPath"
          }
        },
        {
          "name": "compartment sizes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment",
            "grammar": "XPath"
          }
        }
      ],
      "availableSoftwareInterfaceTypes": [
        "library",
        "command-line application",
        "desktop application",
        "BioSimulators Docker image"
      ],
      "dependencies": [
        {
          "name": "libRoadRunner",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "http://libroadrunner.org/"
        },
        {
          "name": "SUNDIALS",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "https://computing.llnl.gov/projects/sundials"
        }
      ]
    },
    {
      "id": "euler",
      "name": "Forward Euler method",
      "kisaoId": {
        "namespace": "KISAO",
        "id": "KISAO_0000030"
      },
      "modelingFrameworks": [
        {
          "namespace": "SBO",
          "id": "SBO_0000293"
        }
      ],
      "modelFormats": [
        {
          "namespace": "EDAM",
          "id": "format_2585",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "modelChangePatterns": [
        {
          "name": "Change component attributes",
          "types": [
            "SedAttributeModelChange",
            "SedComputeAttributeChangeModelChange",
            "SedSetValueAttributeModelChange"
          ],
          "target": {
            "value": "//*/@*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Add components",
          "types": [
            "SedAddXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Remove components",
          "types": [
            "SedRemoveXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Change components",
          "types": [
            "SedChangeXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        }
      ],
      "simulationFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3685",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "simulationTypes": [
        "SedUniformTimeCourseSimulation"
      ],
      "archiveFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3686",
          "version": "1",
          "supportedFeatures": []
        }
      ],
      "parameters": [],
      "outputDimensions": [
        {
          "namespace": "SIO",
          "id": "SIO_000418"
        }
      ],
      "outputVariablePatterns": [
        {
          "name": "time",
          "symbol": {
            "value": "time",
            "namespace": "urn:sedml:symbol"
          }
        },
        {
          "name": "species concentrations",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species",
            "grammar": "XPath"
          }
        },
        {
          "name": "parameter values",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter",
            "grammar": "XPath"
          }
        },
        {
          "name": "reaction fluxes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction",
            "grammar": "XPath"
          }
        },
        {
          "name": "compartment sizes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment",
            "grammar": "XPath"
          }
        }
      ],
      "availableSoftwareInterfaceTypes": [
        "library",
        "command-line application",
        "desktop application",
        "BioSimulators Docker image"
      ],
      "dependencies": [
        {
          "name": "libRoadRunner",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "http://libroadrunner.org/"
        }
      ],
      "citations": []
    },
    {
      "id": "runge_kutta_4",
      "name": "Runge-Kutta fourth order method",
      "kisaoId": {
        "namespace": "KISAO",
        "id": "KISAO_0000032"
      },
      "modelingFrameworks": [
        {
          "namespace": "SBO",
          "id": "SBO_0000293"
        }
      ],
      "modelFormats": [
        {
          "namespace": "EDAM",
          "id": "format_2585",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "modelChangePatterns": [
        {
          "name": "Change component attributes",
          "types": [
            "SedAttributeModelChange",
            "SedComputeAttributeChangeModelChange",
            "SedSetValueAttributeModelChange"
          ],
          "target": {
            "value": "//*/@*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Add components",
          "types": [
            "SedAddXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Remove components",
          "types": [
            "SedRemoveXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Change components",
          "types": [
            "SedChangeXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        }
      ],
      "simulationFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3685",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "simulationTypes": [
        "SedUniformTimeCourseSimulation"
      ],
      "archiveFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3686",
          "version": "1",
          "supportedFeatures": []
        }
      ],
      "parameters": [],
      "outputDimensions": [
        {
          "namespace": "SIO",
          "id": "SIO_000418"
        }
      ],
      "outputVariablePatterns": [
        {
          "name": "time",
          "symbol": {
            "value": "time",
            "namespace": "urn:sedml:symbol"
          }
        },
        {
          "name": "species concentrations",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species",
            "grammar": "XPath"
          }
        },
        {
          "name": "parameter values",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter",
            "grammar": "XPath"
          }
        },
        {
          "name": "reaction fluxes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction",
            "grammar": "XPath"
          }
        },
        {
          "name": "compartment sizes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment",
            "grammar": "XPath"
          }
        }
      ],
      "availableSoftwareInterfaceTypes": [
        "library",
        "command-line application",
        "desktop application",
        "BioSimulators Docker image"
      ],
      "dependencies": [
        {
          "name": "libRoadRunner",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "http://libroadrunner.org/"
        }
      ],
      "citations": []
    },
    {
      "id": "runge_kutta_45",
      "name": "Fehlberg method",
      "kisaoId": {
        "namespace": "KISAO",
        "id": "KISAO_0000086"
      },
      "modelingFrameworks": [
        {
          "namespace": "SBO",
          "id": "SBO_0000293"
        }
      ],
      "modelFormats": [
        {
          "namespace": "EDAM",
          "id": "format_2585",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "modelChangePatterns": [
        {
          "name": "Change component attributes",
          "types": [
            "SedAttributeModelChange",
            "SedComputeAttributeChangeModelChange",
            "SedSetValueAttributeModelChange"
          ],
          "target": {
            "value": "//*/@*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Add components",
          "types": [
            "SedAddXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Remove components",
          "types": [
            "SedRemoveXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Change components",
          "types": [
            "SedChangeXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        }
      ],
      "simulationFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3685",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "simulationTypes": [
        "SedUniformTimeCourseSimulation"
      ],
      "archiveFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3686",
          "version": "1",
          "supportedFeatures": []
        }
      ],
      "parameters": [
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000485"
          },
          "id": "minimum_time_step",
          "name": "Minimum time step",
          "type": "float",
          "value": "1e-12",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000467"
          },
          "id": "maximum_time_step",
          "name": "Maximum time step",
          "type": "float",
          "value": "1.0",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000597"
          },
          "id": "epsilon",
          "name": "Epsilon",
          "type": "float",
          "value": "0.000000000001",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        }
      ],
      "outputDimensions": [
        {
          "namespace": "SIO",
          "id": "SIO_000418"
        }
      ],
      "outputVariablePatterns": [
        {
          "name": "time",
          "symbol": {
            "value": "time",
            "namespace": "urn:sedml:symbol"
          }
        },
        {
          "name": "species concentrations",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species",
            "grammar": "XPath"
          }
        },
        {
          "name": "parameter values",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter",
            "grammar": "XPath"
          }
        },
        {
          "name": "reaction fluxes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction",
            "grammar": "XPath"
          }
        },
        {
          "name": "compartment sizes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment",
            "grammar": "XPath"
          }
        }
      ],
      "availableSoftwareInterfaceTypes": [
        "library",
        "command-line application",
        "desktop application",
        "BioSimulators Docker image"
      ],
      "dependencies": [
        {
          "name": "libRoadRunner",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "http://libroadrunner.org/"
        }
      ],
      "citations": []
    },
    {
      "id": "gillespie",
      "name": "Gillespie direct method of the Stochastic Simulation Algorithm (SSA)",
      "kisaoId": {
        "namespace": "KISAO",
        "id": "KISAO_0000029"
      },
      "modelingFrameworks": [
        {
          "namespace": "SBO",
          "id": "SBO_0000295"
        }
      ],
      "modelFormats": [
        {
          "namespace": "EDAM",
          "id": "format_2585",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "modelChangePatterns": [
        {
          "name": "Change component attributes",
          "types": [
            "SedAttributeModelChange",
            "SedComputeAttributeChangeModelChange",
            "SedSetValueAttributeModelChange"
          ],
          "target": {
            "value": "//*/@*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Add components",
          "types": [
            "SedAddXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Remove components",
          "types": [
            "SedRemoveXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Change components",
          "types": [
            "SedChangeXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        }
      ],
      "simulationFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3685",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "simulationTypes": [
        "SedUniformTimeCourseSimulation"
      ],
      "archiveFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3686",
          "version": "1",
          "supportedFeatures": []
        }
      ],
      "citations": [
        {
          "title": "Exact stochastic simulation of coupled chemical reactions",
          "authors": "Daniel T. Gillespie",
          "journal": "Journal of Physical Chemistry",
          "volume": "81",
          "issue": "25",
          "pages": "2340-2361",
          "year": 1977,
          "identifiers": [
            {
              "namespace": "doi",
              "id": "10.1021/j100540a008",
              "url": "https://doi.org/10.1021/j100540a008"
            }
          ]
        }
      ],
      "parameters": [
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000488"
          },
          "id": "seed",
          "name": "Random number generator seed",
          "type": "integer",
          "value": null,
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000673"
          },
          "id": "nonnegative",
          "name": "Skip reactions which would result in negative species amounts",
          "type": "boolean",
          "value": "false",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        }
      ],
      "outputDimensions": [
        {
          "namespace": "SIO",
          "id": "SIO_000418"
        }
      ],
      "outputVariablePatterns": [
        {
          "name": "time",
          "symbol": {
            "value": "time",
            "namespace": "urn:sedml:symbol"
          }
        },
        {
          "name": "species counts",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species",
            "grammar": "XPath"
          }
        },
        {
          "name": "parameter values",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter",
            "grammar": "XPath"
          }
        },
        {
          "name": "reaction fluxes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction",
            "grammar": "XPath"
          }
        },
        {
          "name": "compartment sizes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment",
            "grammar": "XPath"
          }
        }
      ],
      "availableSoftwareInterfaceTypes": [
        "library",
        "command-line application",
        "desktop application",
        "BioSimulators Docker image"
      ],
      "dependencies": [
        {
          "name": "libRoadRunner",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "http://libroadrunner.org/"
        }
      ]
    },
    {
      "id": "nleq2",
      "name": "Newton-type method for solving non-linear (NL) equations (EQ)",
      "kisaoId": {
        "namespace": "KISAO",
        "id": "KISAO_0000569"
      },
      "modelingFrameworks": [
        {
          "namespace": "SBO",
          "id": "SBO_0000293"
        }
      ],
      "modelFormats": [
        {
          "namespace": "EDAM",
          "id": "format_2585",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "modelChangePatterns": [
        {
          "name": "Change component attributes",
          "types": [
            "SedAttributeModelChange",
            "SedComputeAttributeChangeModelChange",
            "SedSetValueAttributeModelChange"
          ],
          "target": {
            "value": "//*/@*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Add components",
          "types": [
            "SedAddXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Remove components",
          "types": [
            "SedRemoveXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        },
        {
          "name": "Change components",
          "types": [
            "SedChangeXmlModelChange"
          ],
          "target": {
            "value": "//*",
            "grammar": "XPath"
          }
        }
      ],
      "simulationFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3685",
          "version": null,
          "supportedFeatures": []
        }
      ],
      "simulationTypes": [
        "SedSteadyStateSimulation"
      ],
      "archiveFormats": [
        {
          "namespace": "EDAM",
          "id": "format_3686",
          "version": "1",
          "supportedFeatures": []
        }
      ],
      "citations": [
        {
          "title": "A family of Newton codes for systems of highly nonlinear equations - algorithm, implementation, application",
          "authors": "Ulrich Nowak & Lutz Weimann",
          "journal": "Konrad-Zuse-Zentrum für Informationstechnik Berlin",
          "volume": "91-10",
          "year": 1991,
          "identifiers": [
            {
              "namespace": "citeseerx",
              "id": "10.1.1.43.3751",
              "url": "http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.43.3751"
            }
          ]
        },
        {
          "title": "Newton methods for nonlinear problems",
          "authors": "Peter Deuflhard",
          "journal": "Affine Invariance and Adaptive Algorithms",
          "year": 2004,
          "identifiers": [
            {
              "namespace": "doi",
              "id": "10.1007/978-3-642-23899-4",
              "url": "https://doi.org/10.1007/978-3-642-23899-4"
            }
          ]
        }
      ],
      "parameters": [
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000209"
          },
          "id": "relative_tolerance",
          "name": "Relative tolerance",
          "type": "float",
          "value": "1e-12",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000486"
          },
          "id": "maximum_iterations",
          "name": "Maximum number of iterations",
          "type": "integer",
          "value": "100",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000487"
          },
          "id": "minimum_damping",
          "name": "Minimum damping factor",
          "type": "float",
          "value": "1e-20",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000674"
          },
          "id": "allow_presimulation",
          "name": "Presimulate",
          "type": "boolean",
          "value": "false",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000677"
          },
          "id": "presimulation_maximum_steps",
          "name": "Maximum number of steps for presimulation",
          "type": "integer",
          "value": "100",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000680"
          },
          "id": "presimulation_time",
          "name": "Amount of time to pre-simulate",
          "type": "float",
          "value": "100",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000682"
          },
          "id": "allow_approx",
          "name": "Whether to find an approximate solution if an exact solution could not be found",
          "type": "boolean",
          "value": "false",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000683"
          },
          "id": "approx_tolerance",
          "name": "Tolerance for finding an approximate solution",
          "type": "float",
          "value": "0.000001",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000678"
          },
          "id": "approx_maximum_steps",
          "name": "Maximum number of steps for finding an approximate solution",
          "type": "integer",
          "value": "10000",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000679"
          },
          "id": "approx_time",
          "name": "Maximum amount of time to spend finding finding an approximate solution",
          "type": "float",
          "value": "10000",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000675"
          },
          "id": "broyden_method",
          "name": "Broyden method",
          "type": "integer",
          "value": "0",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        },
        {
          "kisaoId": {
            "namespace": "KISAO",
            "id": "KISAO_0000676"
          },
          "id": "linearity",
          "name": "Degree of linearity of the system",
          "type": "integer",
          "value": "3",
          "recommendedRange": null,
          "availableSoftwareInterfaceTypes": [
            "library",
            "command-line application",
            "desktop application",
            "BioSimulators Docker image"
          ]
        }
      ],
      "outputDimensions": [],
      "outputVariablePatterns": [
        {
          "name": "species concentrations",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species",
            "grammar": "XPath"
          }
        },
        {
          "name": "parameter values",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter",
            "grammar": "XPath"
          }
        },
        {
          "name": "reaction fluxes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfReactions/sbml:reaction",
            "grammar": "XPath"
          }
        },
        {
          "name": "compartment sizes",
          "target": {
            "value": "/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment",
            "grammar": "XPath"
          }
        }
      ],
      "availableSoftwareInterfaceTypes": [
        "library",
        "command-line application",
        "desktop application",
        "BioSimulators Docker image"
      ],
      "dependencies": [
        {
          "name": "libRoadRunner",
          "version": null,
          "required": true,
          "freeNonCommercialLicense": true,
          "url": "http://libroadrunner.org/"
        }
      ]
    }
  ],
  "interfaceTypes": [
    "library",
    "command-line application",
    "desktop application",
    "BioSimulators Docker image"
  ],
  "supportedOperatingSystemTypes": [
    "platform-independent"
  ],
  "supportedProgrammingLanguages": [
    {
      "namespace": "Linguist",
      "id": "Python"
    },
    {
      "namespace": "Linguist",
      "id": "C"
    },
    {
      "namespace": "Linguist",
      "id": "C++"
    },
    {
      "namespace": "Linguist",
      "id": "C#"
    }
  ],
  "funding": [
    {
      "funder": {
        "namespace": "FunderRegistry",
        "id": "100000057"
      },
      "grant": "R01-GM123032",
      "url": "https://grantome.com/grant/NIH/R01-GM123032"
    },
    {
      "funder": {
        "namespace": "FunderRegistry",
        "id": "100000057"
      },
      "grant": "R01-GM081070",
      "url": "https://grantome.com/grant/NIH/R01-GM081070"
    }
  ],
  "biosimulators": {
    "specificationVersion": "1.0.0",
    "imageVersion": "1.0.0",
    "validated": false,
    "validationTests": null
  }
}
```
# Guidelines for using Vega to visualize simulation results


## Overview

We recommend [Vega](https://vega.github.io/vega/) for data visualizations of simulation results. Vega is a powerful, declarative grammar for describing interactive, two-dimensional data visualizations.

One key feature of Vega is that it modularly captures the graphical marks which comprise visualizations and how those graphical marks should be painted with data. This feature makes it easy to produce data visualizations for multiple simulation conditions by reusing the same graphical marks with results from multiple simulations. This feature also makes the provenance of data visualizations transparent. As a result, we believe Vega is ideal for collaboration and publication.

Below, we provide recommendations for using Vega to visualize the results of simulation experiments described with SED-ML.

## Creating Vega visualizations of the results of SED-ML files in COMBINE archives

BioSimulators recommends using Vega visualizations with SED-ML as follows:

1. Annotate the Vega signals whose values should be rendered with the values of attributes of simulations or reports of SED-ML documents (e.g., number of a steps of a uniform time course simulation).

    - To set the `value` attribute of a Vega signal equal to the value of an attribute of a simulation or report of a SED-ML document, add a `sedmlUri` key to the signal with a value equal to a list of the location of the SED-ML document, the id of the SED-ML simulation or report, and the name of the attribute of the simulation or report (e.g., `['location/of/simulation.sedml', 'simulationId', 'numberOfSteps']`). To indicate that a signal should be rendered with a list of the values of an attribute of multiple simulations or reports, use `SedDocument:*`, `Simulation:*`, or `Report:*` for the SED-ML document location or simulation/report id (e.g., `['SedDocument:*', 'Report:*', 'id']` to render a signal with a list of the ids of the all of the reports of all of the SED-ML files in the parent COMBINE/OMEX archive).
    - Similarly, to set the `bind` attribute of a Vega signal equal to the value of an attribute of a simulation or report of a SED-ML document, add a `sedmlUri` key to the `bind` attribute with a value as described above.

2. Annotate the Vega data sets whose values should be rendered with the results of SED-ML reports by adding `sedmlUri` keys to these Vega data sets. The values of these keys should be set as follows to indicate the simulation results that should be linked to each Vega data set:
    - To render a Vega data set with the results of all reports from all of the SED-ML files in the parent COMBINE/OMEX archive, the value of the `sedmlUri` key should be an empty array (i.e. `[]`).
    - To render a Vega data set with the result of a single report from one SED-ML file in the parent COMBINE/OMEX archive, the value of the `sedmlUri` key should be a list of the location of the SED-ML document and the id of the report in the document (e.g., `['location/of/simulation.sedml', 'reportId']`).

3. Use the URI `http://purl.org/NET/mediatypes/application/vnd.vega.v5+json` to indicate the formats of Vega files in the manifests of COMBINE/OMEX archives.


### Example data visualization snippet (Vega document that indicates which Vega data sets should be mapped to SED-ML reports)
```json
{
  "signals": [
    {
      "name": "Regular Vega signal",
      "value": 2
    },
    {
      "name": "Vega signal whose value should be rendered with the value of an attribute of a simulation of a SED-ML document",
      "sedmlUri": [
        "simulation_1.sedml",
        "simulation_1",
        "numberOfSteps"
      ]
    }
  ],
  "data": [
    {
      "name": "Regular Vega data set, such as for data for visual marks",
      "values": []
    },
    {
      "name": "Vega data set whose value should be rendered with the result of a report of a SED-ML document",
      "sedmlUri": [
        "simulation.sedml",
        "reportId"
      ]
    }
  ]
}
```

## Rendering Vega visualizations of the results of SED-ML files in COMBINE archives


Simulation software tools should render such Vega visualizations linked to SED-ML files in COMBINE/OMEX archives as follows:

1. Execute the SED-ML files in the COMBINE/OMEX archive and save each report.

2. Use the manifest of the archive to identify the Vega visualizations in the archive (contents with the format `http://purl.org/NET/mediatypes/application/vnd.vega.v5+json`).

3. Extract these Vega visualization files from the archive.

4. Use a JSON library to parse the Vega visualization files.

5. Identify the Vega signals whose values should be rendered with the values of attributes of simulations or reports in SED-ML documents (i.e. Vega signals that have `sedmlUri` keys).

6. Set the values of the Vega signals identified in the previous step to the indicated values of attributes of SED-ML simulations or reports. As illustrated below, SED-ML reports should be encoded as lists of objects that represent the results of each SED-ML dataset. Data sets with multidimensional values should be captured using nested lists.

7. Identify the Vega data sets whose values should be rendered with the values of the results of SED-ML reports (i.e. Vega data sets that have `sedmlUri` keys).

8. Set the values of the Vega data sets identified in the previous step to the results of the indicated SED-ML reports. As illustrated below, SED-ML reports should be encoded as lists of objects that represent the results of each SED-ML dataset. Data sets with multidimensional values should be captured using nested lists.

9. Use [Vega-Embed](https://github.com/vega/vega-embed) to render the resultant Vega visualizations.

### Example simulation results (SED-ML report)

```json
[
  {
    "id": "data_gen_mass",
    "label": "data_gen_mass",
    "name": "mass",
    "values": [
      0.80224854,
      0.8050613294000692,
      0.8078839808114837,
      0.8107165288117656,
      0.8135590080996636,
      0.816411453495585
    ]
  },
  {
    "id": "data_gen_time",
    "label": "data_gen_time",
    "name": "time",
    "values": [
      0,
      0.7,
      1.4,
      2.0999999999999996,
      2.8,
      3.5
    ]
  }
]
```

### Example rendered Vega document (Vega document with data sets replaced with the values of SED-ML reports)


```json
{
  "signals": [
    {
      "name": "Regular Vega signal",
      "value": 2
    },
    {
      "name": "Vega signal whose value should be rendered with the value of a simulation attribute of a SED-ML document",
      "value": 100
    }
  ],
  "data": [
    {
      "name": "Regular Vega data set, such as for data for visual marks",
      "values": []
    },
    {
      "name": "Vega data set whose value should be rendered with the result of a report of a SED-ML document",
      "values": [
        {
          "id": "data_gen_mass",
          "label": "data_gen_mass",
          "name": "mass",
          "values": [
            0.80224854,
            0.8050613294000692,
            0.8078839808114837,
            0.8107165288117656,
            0.8135590080996636,
            0.816411453495585
          ]
        },
        {
          "id": "data_gen_time",
          "label": "data_gen_time",
          "name": "time",
          "values": [
            0,
            0.7,
            1.4,
            2.0999999999999996,
            2.8,
            3.5
          ]
        }
      ]
    }
  ]
}
```

## Example COMBINE/OMEX archives with Vega visualizations

--8<--
vega-examples.md
--8<--

## Tools for converting visualizations of models into Vega data visualizations

--8<--
vega-converters.md
--8<--

## Recommended resources for creating and rendering visualizations

--8<--
vega-resources.md
--8<--# Conventions for describing reproducible and reusable simulation experiments with SED-ML

## Overview

Simulators should support SED-ML L1V3 or later. To accommodate a wide range of modeling frameworks and simulation algorithms, BioSimulators and BioSimulations embrace the additional conventions for SED-ML described below, as well as the conventions for executing SED-ML documents described [here](./simulator-interfaces.md).

## Model and data descriptor source paths

SED-ML can refer to model and data descriptor files in multiple ways, including via paths to local files, URLs, URI fragments to other models defined in the same SED-ML document, and identifiers for an Identifiers.org namespace such as BioModels. When referencing files via local paths, SED-ML documents should use paths relative to the location of the SED-ML document.

To ensure that COMBINE/OMEX archives are self-contained, we encourage SED-ML documents in COMBINE/OMEX archives to reference files via relative paths within archives or other models within the same SED-document.

## Concrete XPath targets for changes to XML-encoded models

SED-ML enables investigators to use [XPaths](https://en.wikipedia.org/wiki/XPath) to specify changes to models that are encoded in XML files. This encompasses models described using CellML, SBML, and other languages. SED-ML documents should use valid XPaths that resolve to XML elements. For example, `/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='A']/@initialConcentration` could be used to indicate a change to the initial condition of the species with id `A`.

In addition, the namespace prefixes used in XPaths should be defined within the SED-ML document as illustrated below.

```xml
<sedML xmlns:sbml="http://www.sbml.org/sbml/level3/version1/core">
  <listOfDataGenerators>
    <dataGenerator>
      <listOfVariables>
        <variable target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='A']/@initialConcentration" />
      </listOfVariables>
    </dataGenerator>
  </listOfDataGenerators>
</sedML>
```

!!!note
    The SED-ML L1V3 and earlier specifications suggest that incomplete XPaths such as `/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='A']` should be used to indicate changes to model elements. We discourage this convention of partial XPaths because these XPaths do not point to the attribute that is intended to be changed. We encourage investigators to use complete XPaths.

!!!note
    SED-ML L1V4 and later documents can use `target` and `symbol` together to reference implicit attributes of model elements, such as fluxes of reactions of flux-balance models.

## Namespaces for `NewXML` elements of changes to XML-encoded models

SED-ML documents can use `sedml:newXML` elements of `sedml:addXML` and `sedml:changeXML` elements to specify objects that should be added to models or replaced in models. SED-ML documents should define the namespace(s) of the content of these `NewXML` elements. For example, a parameter that should be added to a SBML model could be described as `<sbml:parameter xmlns:sbml="http://www.sbml.org/sbml/level3/version1" id="NewParameter" value="10.0" />`.

!!!note 
    The SED-ML specifications suggest that namespaces don't need to be defined for the content of `NewXML` elements. We discourage this convention because XML files which embrace this convention are not consistent with SED-ML's XML schema. We encourage investigators to explicitly define the namespaces involved in the content of `NewXML` elements.

## Data types for model attribute changes and algorithm parameters

SED-ML specifies that the new values of model attribute changes (`sedml:changeAttribute/@sedml:newValue`) and values of algorithm parameters (`sedml:algorithmParameter/@sedml:value`) must be encoded into strings. To ensure that SED-ML files are portable across simulation tools, we define several data types for model attribute changes and algorithm parameters and outlines how each data type should be encoded into strings. The data type of each algorithm parameter should be defined in the specification of each simulation tool.

- `boolean`: Represents Boolean values. Should be encoded into strings as `true`/`false` or `0`/`1`.
- `integer`: Represents integers. Should be encoded in decimal notation (e.g., `1234`).
- `float`: Represents floating point numbers. Should be encoded in decimal (e.g., `1234.567`) or scientific (e.g., `1.234567e3`) notation.
- `string`: Represents strings. Requires no additional encoding.
- `kisaoId`: Represents a KiSAO term. Should be encoding using the id of the term (e.g., `KISAO_0000029`).
- `list`: Represents a list of scalar values. Should be encoding using JSON (e.g., `['a', 'b', 'c']` or `[1, 2, 3]`). For example, the value of the deterministic reactions partition (`KISAO_0000534`) of the Pahle hybrid discrete/continuous Fehlberg method (`KISAO_0000563`) should be a list of the ids of the reactions which should be simulated by the Fehlberg sub-method. Its value should be encoded into SED-ML as `<algorithmParameter kisaoID="KISAO:0000534" value='["ReactionId-1", "ReactionId-1", ...]' />`.
- `object`: Represents key-value pairs. Should be encoding using JSON (e.g., `{a: 1, b: 2}` or `{a: 'x', b: 'y'}`).
- `any`: Represents any other data type. Should be encoding using JSON (e.g., `[{a: 1, b: 2}]`).

Enumerations for the value of an algorithm parameter values can be defined in the specification of a simulator using the `recommendedRange` attribute. This can be combined with any of the above data types.

## Limit use of repeated tasks to the execution of independent simulation runs

In addition to capturing multiple independent simulation runs, `sedml:repeatedTask/@resetModel="False"` provides limited abilities to describe sets of dependent simulation runs, where each run begins from the end state of the previous run. This provides investigators limited abilities to describe meta simulation algorithms.

Simulation tools are encouraged to support a simpler subset of the features of `sedml:repeatedTask` that is sufficient to describe multiple independent simulation runs.

- `sedml:repeatedTask`: Simulation tools should support `resetModel="True"` as described in the SED-ML specifications; the model specifications and initial conditions should be reset. Simulator state such as the states of random number generators should not be reset. When `resetModel="False"`, simulation tools should support limited preservation of the state of simulations between iterations. Simulation tools should accumulate changes to the specifications of the model(s) involved in the task. Simulations tools should not copy the final simulation state from the previous iteration to the initial state of the next iteration.

- Sub-tasks (`sedml:subTask`): Successive subtasks should be executed independently, including when they involve the same model. The final state of the previous sub-task should not be used to set up the initial state for the next sub-task.

- Shape of model variables for the results of repeated tasks: Repeated tasks should produce multi-dimensional results. The first dimension should represent the iterations of the main range of the repeated task. The second dimension should represent the sub-tasks of the repeated task. The results of sub-tasks should be ordered in the same order the sub-tasks were executed (in order of their order attributes). The result of each sub-task should be reshaped to the largest shape of its sibling sub-tasks by padding smaller results with `NaN`. Each nesting of repeated tasks should contribute two additional dimensions for their ranges and sub-tasks. The final dimensions should be the dimensions of the atomic tasks of the repeated task (e.g., time for tasks of uniform time courses).

## Canonical order of execution of tasks

For reproducibility, simulation tools should execute tasks in the order in which they are defined in SED-ML files.

Furthermore, because the order of execution can affect the results of simulations, in general, each task should be executed, including tasks which do not contribute to any output. This is particularly important for simulation tools that implement Monte Carlo algorithms. One exception is tasks whose results are invariant to their order of execution, such as most deterministic simulations. Such tasks can be executed in any order or in parallel.

## Limit use of symbols to variables of data generators

SED-ML uses symbols to reference implicit properties of simulations that are not explicitly defined in the specification of the model for the simulation. The most frequently used symbol for SBML-encoded models is `urn:sedml:symbol:time` for the variable time. Such symbols only have defined values for simulations of models and not for models themselves.

Consequently, symbols should only be used in contexts where simulations are defined. Specifically, symbols should only be used in conjunction with variables of `sedml:dataGenerator` to record predicted values of symbols. Symbols should not be used in conjunction with the variables of `sedml:computeChange`, `sedml:setValue`, or `sedml:functionalRange`. Symbols should also not be used with `sedml:setValue` to set the values of symbols.

## Variable targets for model objects that generate multiple predictions

Some algorithms, such as flux balance analysis (FBA, [`KISAO_0000437`](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000437)) and flux variability analysis (FVA, [`KISAO_0000526`](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000526)) generate multiple predictions for each model object. For example, flux variability analysis predicts minimum and maximum fluxes for each reaction. Targets (`sedml:variable/@sedml:target`) for such predictions should indicate the id of the desired prediction. To ensure portability of SED-ML files between simulation tools, we define the following ids. Please use [GitHub issues](https://github.com/biosimulators/Biosimulators/issues/new/choose) to suggest additional ids for additional predictions of other algorithms.

- FBA ([`KISAO_0000437`](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000437)), parsimonious FBA ([`KISAO_0000528`](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000528)), geometric FBA ([`KISAO_0000527`]https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000527):
    - Objective: `fbc:objective/@fbc:value`
    - Reaction flux: `sbml:reaction/@fbc:flux`
    - Reaction reduced cost: `sbml:reaction/@fbc:reducedCost`
    - Species shadow price: `sbml:species/@fbc:shadowPrice`
- FVA ([KISAO_0000526](https://www.ebi.ac.uk/ols/ontologies/kisao/terms?iri=http%3A%2F%2Fwww.biomodels.net%2Fkisao%2FKISAO%23KISAO_0000526)):
    - Minimum reaction flux: `sbml:reaction/@fbc:minFlux`
    - Maximum reaction flux: `sbml:reaction/@fbc:maxFlux`

## Unique data set labels
To facilitate automated interpretation of simulation results, the data sets within a report should have unique labels (`sedml:dataSet/@sedml:label`). Note, the same label can be used across multiple reports.


## Guides for using SED-ML and the COMBINE archive format with specific model languages

Simulation tools should recognize the URNs and IRIs below to identify model languages described in SED-ML files and COMBINE/OMEX archives. The links in the "Info" column below contain more information about how simulation tools should interpret SED-ML in combination with specific model languages.

| Language                                                              | EDAM id | SED-ML URN	               | COMBINE archive specification URI	                   | MIME type              | Extensions     |Info                                                                           |
| ---                                                                   | ---     | ---                        | ---                                                  |---                     | ---            |---                                                                            |
| [BNGL](https://bionetgen.org/)                                        | 3972    | urn:sedml:language:bngl    | http://purl.org/NET/mediatypes/text/bngl+plain	   | text/bngl+plain	    | .bngl	         | [:link:](https://docs.biosimulators.org/Biosimulators_BioNetGen/tutorial.html)| 
| [CellML](https://cellml.org/)                                         | 3240    | urn:sedml:language:cellml  | http://identifiers.org/combine.specifications/cellml | application/cellml+xml | .xml, .cellml	 | [:link:](https://sed-ml.org/specifications.html)                              | 
| ([NeuroML](https://neuroml.org/))/[LEMS](https://lems.github.io/LEMS/)| 9004    | urn:sedml:language:lems    | http://purl.org/NET/mediatypes/application/lems+xml  | application/lems+xml   | .xml	         |                                                                               |
| [SBML](http://sbml.org/)                                              | 2585    | urn:sedml:language:sbml	   | http://identifiers.org/combine.specifications/sbml   | application/sbml+xml   | .xml, .sbml	 | [:link:](https://sed-ml.org/specifications.html)                              |
| [Smoldyn](http://www.smoldyn.org/)                                    | 9001    | urn:sedml:language:smoldyn | http://purl.org/NET/mediatypes/text/smoldyn+plain	   | text/smoldyn+plain     | .txt	         | [:link:](https://github.com/ssandrews/Smoldyn/blob/master/Using-Smoldyn-with-SED-ML-COMBINE-BioSimulators.md)     |

Example SED-ML files and COMBINE/OMEX archives for all of the languages listed above are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples#compatibility-of-the-example-archives-with-simulation-tools).

## Recommended resources for implementing the execution of simulation experiments

Below are helpful tools for implementing the execution of simulation experiments described with SED-ML:

- [BioSimulators utils ](https://docs.biosimulators.org/Biosimulators_utils/) is a Python library which provides functions for implementing command-line interfaces to the above specifications, as well as functions for interpreting COMBINE/OMEX archives and SED-ML files, generating tables and plots of simulation plots, and logging the execution of COMBINE/OMEX archives. BioSimulators utils provides high-level access to some of the lower-level libraries listed below.

- [libSED-ML](https://github.com/fbergmann/libSEDML) is a library for serializing and deserializing SED-ML documents to and from XML files. libSED-ML provides bindings for several languages.

- [jlibSED-ML](https://sourceforge.net/projects/jlibsedml/) is a Java library for serializing and deserializing SED-ML documents to and from XML files. The library also provides methods for resolving models, working with XPath targets for model elements, applying model changes, orchestrating the execution of tasks, calculating the values of data generators, and logging the execution of simulations. Note, jLibSED-ML support SED-ML <= L1V2 and diverges from some of the conventions described here.# Guidelines for Docker images for biosimulation tools

## Overview
The BioSimulators standard for Docker images of biosimulation software tools outlines the syntax and semantics of the entry points of Docker images of biosimulators. This standard ensures that containerized simulation tools can be consistently executed with the same input arguments (a path to a COMBINE/OMEX archive that defines models and simulations and a path to save the outputs of the simulation experiments defined in the archive) and that simulation tools produce consistent outputs (reports and plots at consistent paths in consistent formats). The standard also specifies how images of biosimulation software tools should use Docker labels to capture metadata about themselves.

## Specifications of entry points 

Simulation tools which implement the BioSimulators Docker image standard should provide an entry point that maps to a command-line interface that implements the [BioSimulators standard](./simulator-interfaces.md):

The entry point should map to a command-line interface for a biosimulator (`ENTRYPOINT ["simulator-standard-interface"]`). The default arguments for the entry point should be an empty list (`CMD []`).

## Specifications of Docker labels

Simulation tools which implement the BioSimulators Docker image standard should use labels to provide metadata consistent with the [Open Containers Initiative](https://opencontainers.org/) and [BioContainers](https://biocontainers.pro/).

Specifically, images should provide the following labels:

- Open Containers Initiative labels:
    - `org.opencontainers.image.title`: Human-readable title of the image.
    - `org.opencontainers.image.version`: Version of the software in the image.
    - `org.opencontainers.image.revision`: Source control revision identifier of the software in the image.
    - `org.opencontainers.image.description`: Human-readable description of the software in the image.
    - `org.opencontainers.image.url`: URL to find more information about the image.
    - `org.opencontainers.image.documentation`: URL to get documentation about the image.
    - `org.opencontainers.image.source`: URL to get the Dockerfile for building the image.
    - `org.opencontainers.image.authors`: Contact details of the people or organization responsible for the image.
    - `org.opencontainers.image.vendor`: Name of the entity, organization or individual which distributes the image.
    - `org.opencontainers.image.licenses`: [SPDX](https://spdx.org/) expression which describes the license(s) under which the software in the image is distributed.
    - `org.opencontainers.image.created`: Date and time when the image was built (RFC 3339).

- BioContainers labels:
    - `version`: Version of the image (e.g., `1.0.0`)
    - `software`: Simulation program wrapped into the image (e.g., `BioNetGen`).
    - `software.version`: Version of the simulation program wrapped into the image (e.g., `2.5.0`).
    - `about.summary`: Short description of the simulation program (e.g., `Package for rule-based modeling of complex biochemical systems`).
    - `about.home`: URL for the simulation program (e.g., `https://bionetgen.org/`).
    - `about.documentation`: URL for documentation for the simulation program (e.g., `https://bionetgen.org/`).
    - `about.license_file`: URL for the license for the simulation program (e.g., `https://github.com/RuleWorld/bionetgen/blob/master/LICENSE`).
    - `about.license`: SPDX license id for the license for the simulation program (e.g., `SPDX:MIT`). See [SPDX](https://spdx.org/) for a list of licenses and their ids.
    - `about.tags`: Comma-separated list of tags which describe the simulation program (e.g., `rule-based modeling`,`dynamical simulation`,`systems biology`,`BNGL`,`BioSimulators`). Please include the tag `BioSimulators`.
    - `extra.identifiers.biotools`: Optionally, the bio.tools identifier for the simulation program (e.g., `bionetgen`). Visit [bio.tools](https://bio.tools/) to request identifiers for simulation programs.
    - `maintainer`: Name and email of the person/team who developed the image (e.g., `Jonathan Karr <karr@mssm.edu>`).'


## Additional recommendations for best practices

To ensure that containerized simulation tools can be executed inside high-performance computing clusters, where root access is typically not allowed and conversion to Singularity images is necessary, we recommend that developers also follow the best practices below for Dockerfiles. For more discussion, we recommend [Syslab's best practice guide](https://sylabs.io/guides/3.7/user-guide/singularity_and_docker.html#best-practices).

- Sources of containerized simulation tools and their dependencies: To ensure that the construction of Docker images is reproducible and portable, the simulation tools inside images should be installed from internet sources rather than the local file system. One exception is licenses that are needed to install commercial software. These can be copied from a local directory such as `assets/`, deleted and squashed out of the final image, and injected again when the image is executed.
- Installation locations of containerized simulation tools and their dependencies: Because Docker images are typically run as root, `/root` should be reserved for the home directory of the user which executes the image. Similarly, `/tmp` should be reserved for temporary files that must be created during the execution of the image. Consequently, the simulation tools inside containers and their dependencies should be installed to different directories other than `/root` and `/tmp`.
- Environment variables: All environment variables that the containerized simulation tool supports should be explicitly defined using the `ENV` directive.
User privileges: Do not use the `USER` directive.

## Example

An example Dockerfile for [tellurium](http://tellurium.analogmachine.org/) is below.

```Dockerfile
# Base OS
FROM python:3.9-slim-buster

# metadata
LABEL \
    org.opencontainers.image.title="tellurium" \
    org.opencontainers.image.version="2.1.5" \
    org.opencontainers.image.description="Python-based environment for model building, simulation, and analysis that facilitates reproducibility of models in systems and synthetic biology" \
    org.opencontainers.image.url="http://tellurium.analogmachine.org/" \
    org.opencontainers.image.documentation="https://tellurium.readthedocs.io/" \
    org.opencontainers.image.source="https://github.com/biosimulators/Biosimulators_tellurium" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="Apache-2.0" \
    \
    base_image="python:3.9-slim-buster" \
    version="0.0.1" \
    software="tellurium" \
    software.version="2.1.6" \
    about.summary="Python-based environment for model building, simulation, and analysis that facilitates reproducibility of models in systems and synthetic biology" \
    about.home="http://tellurium.analogmachine.org/" \
    about.documentation="https://tellurium.readthedocs.io/" \
    about.license_file="https://github.com/sys-bio/tellurium/blob/develop/LICENSE.txt" \
    about.license="SPDX:Apache-2.0" \
    about.tags="kinetic modeling,dynamical simulation,systems biology,biochemical networks,SBML,SED-ML,COMBINE,OMEX,BioSimulators" \
    maintainer="BioSimulators Team <info@biosimulators.org>"

# Install requirements
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        libxml2 \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_tellurium
RUN pip install /root/Biosimulators_tellurium

# Entrypoint
ENTRYPOINT ["tellurium"]
CMD []
```
# Format for metadata about simulation projects

## Overview

Metadata about COMBINE/OMEX archives should be annotated using triples of subjects, predicates, and objects in RDF XML files according to the [OMEX Metadata guidelines](https://doi.org/10.1515/jib-2021-0020).

On top of these guidelines, we recommend the predicates and identifier namespaces described below. In addition, the URI for each object should be annotated using `http://dublincore.org/specifications/dublin-core/dcmi-terms/identifier` and a human-readable description of each object should be annotated using `http://www.w3.org/2000/01/rdf-schema#label`. Futhermore, alternative predicates should also be described using `http://dublincore.org/specifications/dublin-core/dcmi-terms/description`.

## Recommended URIs for COMBINE/OMEX archives and the contents

COMBINE/OMEX archives should be referenced using unique identifiers which begin with the prefix `http://omex-library.org/` and end with the file extension `.omex` (e.g., `http://omex-library.org/BioSimulations-0001.omex`). We recommend using identifiers that are a concatenation of `http://omex-library.org/` and the local filename of the COMBINE archive (e.g., `BioSimulations-0001.omex`).

Files in COMBINE archives should be referenced by concatenating the above identifiers with their location within their parent COMBINE archives (e.g., `http://omex-library.org/BioSimulations-0001.omex/simulation.sedml`).

Elements in SED-ML files in COMBINE archives should be referenced by concatenating the above identifiers with their id (e.g., `http://omex-library.org/BioSimulations-0001.omex/simulation.sedml/Figure_3a`).

## Recommended predicates and objects for annotating COMBINE archives

We recommend that COMBINE/OMEX archives and components of archives be annotated using the predicates and objects outlined below.

- Title:
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/title`
    - Object: Literal string
- Abstract (short summary):
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/abstract`
    - Object: Literal string
- Description (long summary):
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/description`
    - Object: String formatted using [GitHub-flavored markdown](https://github.github.com/gfm/) (GFM)
    !!! warning
        BioSimulations uses [Marked](https://marked.js.org/) to render GFM. Marked supports most, but not all features of GFM. More information about Marked's support for GFM is available [here](https://github.com/markedjs/marked/discussions/1202).
- Keyword:
    - Predicate: `http://prismstandard.org/namespaces/basic/2.0/keyword`
    - Object: Literal string
- Thumbnail image:
    - Predicate: `http://www.collex.org/schema#thumbnail`
    - Object: URI for a GIF, JPEG, PNG, or WEBP file inside the COMBINE/OMEX archive (e.g., `http://omex-library.org/Ciliberto.omex/BioSim0001.png`)
        - Format: GIF, JPEG, PNG, or WEBP
        - Size: At least 1216 pixels wide, and readable at approximately 350 pixels wide. Thumbnails are displayed at approximately 352-552 pixels in the project browse view. Thumbnails are displayed at approximately 352-1216 pixels in the project views.
        - Aspect ratio: The optimal aspect ratio  for the project browse view is 1.625
        - File size: No limit, images are automatically optimized
- Organism captured by a project or component of a project:
    - Predicate: `http://biomodels.net/biology-qualifiers/hasTaxon`
    - Objects: Identifiers.org URI for an entry in NCBI Taxonomy (e.g., `http://identifiers.org/taxonomy/9606`), Literal string
- Other biology (e.g., cell type, organ) captured by a project or component of a project:
    - Predicate: `http://biomodels.net/biology-qualifiers/encodes`
    - Objects: URI (e.g., `https://www.uniprot.org/uniprot/P07527`), Literal string
- Source of a project or component of a project (e.g., GitHub repository):
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/source`
    - Objects: URI (e.g., `https://github.com/org/repo`), Literal string
- Predecessor of a project or component of a project:
    - Predicate: `http://biomodels.net/model-qualifiers/isDerivedFrom`
    - Objects: URI (e.g., `http://identifiers.org/biomodels.db:BIOMD0000000296`), Literal string
- Successor of a project or component of a project:
    - Predicate: `http://purl.org/spar/scoro/successor`
    - Objects: URI (e.g., `http://identifiers.org/biomodels.db:BIOMD0000000298`), Literal string
- More information about a project or component of a project:
    - Predicate: `http://www.w3.org/2000/01/rdf-schema#seeAlso`
    - Objects: URI (e.g., `http://mpf.biol.vt.edu/lab_website/`), Literal string
- References for a project or component of a project:
    - Predicate: `http://purl.org/dc/terms/references`
    - Objects: URI (e.g., `http://identifiers.org/pubmed:1234`), Literal string
- Other identifier for a project or component of a project (e.g., in a primary model repository):
    - Predicate: `http://biomodels.net/model-qualifiers/is`
    - Objects: URI (e.g., `http://identifiers.org/biomodels.db:BIOMD0000000297`), Literal string
- Citation for a project or component of a project:
    - Predicate: `http://biomodels.net/model-qualifiers/isDescribedBy`
    - Objects: Identifiers.org DOI URI (e.g., `http://identifiers.org/doi:10.1083/jcb.200306139`, `http://identifiers.org/pubmed:1234`, `http://identifiers.org/arxiv:0807.4956v1`), Literal string
- Author:
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/creator`
    - Objects: ORCID Identifiers.org URI (e.g., `http://identifiers.org/orcid:0000-0001-7560-6013`), Literal string
- Contributor (e.g., curator):
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/contributor`
    - Objects: ORCID Identifiers.org URI (e.g., `http://identifiers.org/orcid:0000-0001-7560-6013`), Literal string
- License for a project or component of a project:
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/identifier`
    - Objects: SPDX Identifiers.org URI (e.g., `http://identifiers.org/spdx:CCO`), Literal string
- Funder:
    - Predicate: `http://purl.org/spar/scoro/funder`
    - Objects: FunderRegistry Identifiers.org URI (e.g., `http://identifiers.org/doi:10.13039/100000185`), Literal string
- Creation date:
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/created`
    - Objects: WC3DTF-encoded literal string (e.g., `2021-06-01`)
- Modification date:
    - Predicate: `http://dublincore.org/specifications/dublin-core/dcmi-terms/modified`
    - Objects: WC3DTF-encoded literal string (e.g., `2021-06-01`)

We also provides the following recommendations:

- The title, abstract, description, license, and creation date predicates should only be used once per subject.

- Identifiers.org URIs should be used when possible because Identifiers.org can be used to validate URIs and map URIs to URLs. For example, `http://identifiers.org/orcid:0000-0001-7560-6013` should be used rather than `https://orcid.org/0000-0001-7560-6013`.

## Required metadata about COMBINE archives for BioSimulations

Submissions to BioSimulations must include the following metadata:

- Title (`http://dublincore.org/specifications/dublin-core/dcmi-terms/title`)

This requirement is currently set low to accommodate old projects in community repositories that have minimal structured metadata. Over time, we aim to raise this requirement.

## Recommendations for describing the SED-ML files and plots responsible for figures
We recommend using the `http://dublincore.org/specifications/dublin-core/dcmi-terms/identifier` predicate and literal strings to describe the SED-ML files, reports, and plots responsible for tables and figures in articles.

```xml
<rdf:Description rdf:about="http://omex-library.org/BioSim0001.omex/sim.sedml/figure1">
  <bqmodel:is>
    <rdf:Description>
      <dc:identifier rdf:resource="https://doi.org/10.1371/journal.pcbi.1008379.g001"/>
      <rdfs:label>Figure 1a</rdfs:label>
    </rdf:Description>
  </bqmodel:is>
</rdf:Description>
```

## Recommendations for describing the provenance of computationally-generated files
We recommend using the `http://biomodels.net/model-qualifiers/isDerivedFrom` predicate to indicate the source of computationally-generated files, such as SED-ML files automatically created from model files (e.g., SBML). The subjects and objects of such triples should be described using OMEX library URIs (e.g., `http://omex-library.org/BioSim0001.omex/simulation.sedml`, `http://omex-library.org/BioSim0001.omex/model.xml`) that represent their location within their parent COMBINE/OMEX archive.

```xml
<rdf:Description rdf:about="http://omex-library.org/BioSim0001.omex/simulation.sedml">
  <bqmodel:isDerivedFrom>
    <rdf:Description>
      <dc:identifier rdf:resource="http://omex-library.org/BioSim0001.omex/model.xml"/>
      <rdfs:label>model</rdfs:label>
    </rdf:Description>
  </bqmodel:isDerivedFrom>
</rdf:Description>
```

## Example metadata about a simulation project (COMBINE/OMEX archive and SED-ML file)

As an example, below is a representation of metadata for the [Ciliberto 2003 model of the budding yeast cell cycle](https://identifiers.org/doi:10.1083/jcb.200306139).

```xml
<rdf:RDF>
  <!-- metadata about a COMBINE/OMEX archive -->
  <rdf:Description rdf:about="http://omex-library.org/BioSim0001.omex">
    <!-- keywords -->
    <prism:keyword>morphogenesis checkpoint</prism:keyword>
    <prism:keyword>G2</prism:keyword>

    <!-- thumbnail image -->
    <collex:thumbnail rdf:resource="http://omex-library.org/BioSim0001.omex/Figure1.png"/>

    <!-- long description -->
    <dc:description>Based on published observations of budding yeast and analogous control signals in fission yeast. The simulations accurately reproduce the phenotypes of a dozen checkpoint mutants. Among other predictions, the model attributes a new role to Hsl1, a kinase known to play a role in Swe1 degradation: Hsl1 must also be indirectly responsible for potent inhibition of Swe1 activity. The model supports the idea that the morphogenesis checkpoint, like other checkpoints, raises the cell size threshold for progression from one phase of the cell cycle to the next.</dc:description>

    <!-- taxon -->
    <bqbiol:hasTaxon>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/taxonomy/4896"/>
        <rdfs:label>Schizosaccharomyces pombe</rdfs:label>
      </rdf:Description>
    </bqbiol:hasTaxon>

    <!-- gene, RNA, protein, cell type, tissue, organ, disease -->
    <bqbiol:encodes>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/GO:0009653"/>
        <rdfs:label>anatomical structure morphogenesis</rdfs:label>
      </rdf:Description>
    </bqbiol:encodes>

    <bqbiol:encodes>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/kegg:ko04111"/>
        <rdfs:label>Cell cycle - yeast</rdfs:label>
      </rdf:Description>
    </bqbiol:encodes>

    <!-- more information (i.e., domain specific information beyond recommended relationships and identifiers) -->
    <swo:SWO_0000001>
      <rdf:Description>
        <dc:identifier rdf:resource="http://www.math.pitt.edu/~bard/xpp/xpp.html"/>
        <rdfs:label>XPP</rdfs:label>
        <dc:description>Software</dc:description>
      </rdf:Description>
    </swo:SWO_0000001>

    <!-- related things to provide hyperlinks to -->
    <rdfs:seeAlso>
      <rdf:Description>
        <dc:identifier rdf:resource="https://www.bioch.ox.ac.uk/research/novak"/>
        <rdfs:label>Novak Lab</rdfs:label>
      </rdf:Description>
    </rdfs:seeAlso>

    <rdfs:seeAlso>
      <rdf:Description>
        <dc:identifier rdf:resource="http://mpf.biol.vt.edu/lab_website/"/>
        <rdfs:label>Tyson Lab</rdfs:label>
      </rdf:Description>
    </rdfs:seeAlso>

    <!-- other identifiers -->
    <bqmodel:is>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/biomodels.db:BIOMD0000000297"/>
        <rdfs:label>biomodels.db:BIOMD0000000297</rdfs:label>
      </rdf:Description>
    </bqmodel:is>

    <bqmodel:is>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/jws:ciliberto"/>
        <rdfs:label>jws:ciliberto</rdfs:label>
      </rdf:Description>
    </bqmodel:is>

    <!-- authors(s) -->
    <dc:creator>
      <rdf:Description>
        <foaf:name>Andrea Ciliberto</foaf:name>
        <rdfs:label>Andrea Ciliberto</rdfs:label>
      </rdf:Description>
    </dc:creator>

    <dc:creator>
      <rdf:Description>
        <foaf:accountName rdf:resource="https://orcid.org/0000-0002-6961-1366"/>
        <foaf:name>Bela Novak</foaf:name>

        <dc:identifier rdf:resource="http://identifiers.org/orcid:0000-0002-6961-1366"/>
        <rdfs:label>Bela Novak</rdfs:label>
      </rdf:Description>
    </dc:creator>

    <dc:creator>
      <rdf:Description>
        <foaf:accountName rdf:resource="https://orcid.org/0000-0001-7560-6013"/>
        <foaf:name>John J. Tyson</foaf:name>

        <dc:identifier rdf:resource="http://identifiers.org/orcid:0000-0001-7560-6013"/>
        <rdfs:label>John J. Tyson</rdfs:label>
      </rdf:Description>
    </dc:creator>

    <!-- contributor(s) -->
    <dc:contributor>
      <rdf:Description>
        <foaf:accountName rdf:resource="https://orcid.org/0000-0002-2605-5080"/>
        <foaf:name>Jonathan Karr</foaf:name>

        <dc:identifier rdf:resource="http://identifiers.org/orcid:0000-0002-2605-5080"/>
        <rdfs:label>Jonathan Karr</rdfs:label>
      </rdf:Description>
    </dc:contributor>

    <!-- citations -->
    <bqmodel:isDescribedBy>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/doi:10.1083/jcb.200306139"/>
        <rdfs:label>Andrea Ciliberto, Bela Novak &amp; John J. Tyson. Mathematical model of the morphogenesis checkpoint in budding yeast. Journal of Cell Biology 163, 6 (2003): 1243-1254.</rdfs:label>
      </rdf:Description>
    </bqmodel:isDescribedBy>

    <!-- license -->
    <dc:license>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/spdx:CC0-1.0"/>
        <rdfs:label>CC0-1.0</rdfs:label>
      </rdf:Description>
    </dc:license>

    <!-- funder -->
    <scoro:funder>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/doi:10.13039/100000185"/>
        <rdfs:label>Defense Advanced Research Projects Agency</rdfs:label>
      </rdf:Description>
    </scoro:funder>

    <scoro:funder>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/doi:10.13039/100000913"/>
        <rdfs:label>James S. McDonnell Foundation</rdfs:label>
      </rdf:Description>
    </scoro:funder>

    <scoro:funder>
      <rdf:Description>
        <dc:identifier rdf:resource="http://identifiers.org/doi:10.13039/501100010024"/>
        <rdfs:label>Hungarian Science Foundation</rdfs:label>
      </rdf:Description>
    </scoro:funder>

    <!-- created -->
    <dc:created>
      <rdf:Description>
        <dc:W3CDTF>2003-12-22</dc:W3CDTF>
      </rdf:Description>
    </dc:created>

    <!-- modified -->
    <dc:modified>
      <dcterms:W3CDTF>
        <rdf:value>2021-06-25</rdf:value>
      </dcterms:W3CDTF>
    </dc:modified>
  </rdf:Description>


  <!-- metadata about components of the archive -->
  <rdf:Description rdf:about="http://omex-library.org/BioSim0001.omex/simulation_1.sedml">
    <!-- links from SED-ML plots to figures of articles -->
    <bqmodel:is>
      <rdf:Description>
        <dc:identifier rdf:resource="https://doi.org/10.1083/jcb.200306139"/>
        <rdfs:label>Figure 3</rdfs:label>
      </rdf:Description>
    </bqmodel:is>
  </rdf:Description>
</rdf:RDF>
```# Contributing to BioSimulations and BioSimulators

We enthusiastically welcome contributions to BioSimulations and BioSimulators!

## Simulation projects
Simulation projects can be contributed to [BioSimulations](https://biosimulations.org). Please see the [User guide](../../users/publishing-projects) for more information.

## Simulation tools
Simulation tools can be contributed to [BioSimulators](https://biosimulators.org). Please see the [User guide](../../users/publishing-tools) for more information.

## Platform
We welcome contributions by GitHub pull requests to the [BioSimulations Git repository](https://github.com/biosimulations/biosimulations). Please see the [Guide to Contributing](../developers/index.md) for information about how to get started. Please also contact the [core developers](mailto:info@biosimulations.org) to coordinate potential contributions.
# Contributors

Numerous individuals and organizations have contributed to BioSimulations and BioSimulators, including the following:

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/AMICI-dev"><img src="https://avatars.githubusercontent.com/u/68919097?v=4?s=100" width="100px;" alt=""/><br /><sub><b>AMICI</b></sub></a><br /><a href="#tool-AMICI-dev" title="Tools">🔧</a></td>
    <td align="center"><a href="https://fun.bio.keio.ac.jp/"><img src="https://avatars.githubusercontent.com/u/1589676?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Akira Funahashi</b></sub></a><br /><a href="#tool-funasoul" title="Tools">🔧</a></td>
    <td align="center"><a href="https://hellix.com/Alan/"><img src="https://avatars.githubusercontent.com/u/602265?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Alan Garny</b></sub></a><br /><a href="#ideas-agarny" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-agarny" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/ajelenak"><img src="https://avatars.githubusercontent.com/u/7267124?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Aleksandar Jelenak</b></sub></a><br /><a href="#tool-ajelenak" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/ASinanSaglam"><img src="https://avatars.githubusercontent.com/u/11724447?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ali Sinan Saglam</b></sub></a><br /><a href="#data-ASinanSaglam" title="Data">🔣</a></td>
    <td align="center"><a href="https://uni-tuebingen.de/en/127116"><img src="https://avatars.githubusercontent.com/u/1740827?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Andreas Dräger</b></sub></a><br /><a href="#tool-draeger" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/AnkitaxPriya"><img src="https://avatars.githubusercontent.com/u/44089458?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ankita</b></sub></a><br /><a href="#data-AnkitaxPriya" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://ankursinha.in/"><img src="https://avatars.githubusercontent.com/u/102575?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ankur Sinha</b></sub></a><br /><a href="#tool-sanjayankur31" title="Tools">🔧</a></td>
    <td align="center"><a href="https://research.pasteur.fr/en/member/anna-zhukova"><img src="https://avatars.githubusercontent.com/u/10465838?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Anna Zhukova</b></sub></a><br /><a href="#data-annazhukova" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/AnneGoelzer"><img src="https://avatars.githubusercontent.com/u/32333634?v=4?s=100" width="100px;" alt=""/><br /><sub><b>AnneGoelzer</b></sub></a><br /><a href="#data-AnneGoelzer" title="Data">🔣</a></td>
    <td align="center"><a href="https://www.mountsinai.org/profiles/arthur-p-goldberg"><img src="https://avatars.githubusercontent.com/u/33882?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Arthur P Goldberg</b></sub></a><br /><a href="#ideas-artgoldberg" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="http://aurelien.naldi.info/"><img src="https://avatars.githubusercontent.com/u/250984?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Aurélien Naldi</b></sub></a><br /><a href="#data-aurelien-naldi" title="Data">🔣</a> <a href="#tool-aurelien-naldi" title="Tools">🔧</a></td>
    <td align="center"><a href="http://bshaikh.com"><img src="https://avatars.githubusercontent.com/u/32490144?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Bilal Shaikh</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=bilalshaikh42" title="Code">💻</a> <a href="https://github.com/biosimulations/biosimulations/commits?author=bilalshaikh42" title="Documentation">📖</a> <a href="#infra-bilalshaikh42" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
    <td align="center"><a href="https://www.ebi.ac.uk/biomodels"><img src="https://avatars.githubusercontent.com/u/74367888?v=4?s=100" width="100px;" alt=""/><br /><sub><b>BioModels</b></sub></a><br /><a href="#data-EBI-BioModels" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://teusinkbruggemanlab.nl/brett-olivier/"><img src="https://avatars.githubusercontent.com/u/5011985?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brett Olivier</b></sub></a><br /><a href="#tool-bgoli" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/briandrawert"><img src="https://avatars.githubusercontent.com/u/1413538?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brian Drawert</b></sub></a><br /><a href="#tool-briandrawert" title="Tools">🔧</a></td>
    <td align="center"><a href="https://briansimulator.org/"><img src="https://avatars.githubusercontent.com/u/2292949?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brian simulator</b></sub></a><br /><a href="#tool-brian-team" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.copasi.org/"><img src="https://avatars.githubusercontent.com/u/1854399?v=4?s=100" width="100px;" alt=""/><br /><sub><b>COPASI</b></sub></a><br /><a href="#tool-copasi" title="Tools">🔧</a></td>
    <td align="center"><a href="https://reproduciblebiomodels.org"><img src="https://avatars.githubusercontent.com/u/70044163?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Center for Reproducible Biomedical Modeling</b></sub></a><br /><a href="#financial-reproducible-biomedical-modeling" title="Financial">💵</a> <a href="#fundingFinding-reproducible-biomedical-modeling" title="Funding Finding">🔍</a> <a href="#projectManagement-reproducible-biomedical-modeling" title="Project Management">📆</a></td>
    <td align="center"><a href="https://github.com/CiaranWelsh"><img src="https://avatars.githubusercontent.com/u/19502680?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ciaran Welsh</b></sub></a><br /><a href="#tool-CiaranWelsh" title="Tools">🔧</a></td>
    <td align="center"><a href="https://claudine-chaouiya.pedaweb.univ-amu.fr/index.html"><img src="https://avatars.githubusercontent.com/u/40125033?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Claudine Chaouiya</b></sub></a><br /><a href="#data-chaouiya" title="Data">🔣</a> <a href="#tool-chaouiya" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/danv61"><img src="https://avatars.githubusercontent.com/u/29076329?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dan Vasilescu</b></sub></a><br /><a href="#tool-danv61" title="Tools">🔧</a></td>
    <td align="center"><a href="https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/5122/index.html"><img src="https://avatars.githubusercontent.com/u/18048784?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Daniel Weindl</b></sub></a><br /><a href="#tool-dweindl" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/dbrnz"><img src="https://avatars.githubusercontent.com/u/239220?v=4?s=100" width="100px;" alt=""/><br /><sub><b>David Brooks</b></sub></a><br /><a href="#tool-dbrnz" title="Tools">🔧</a></td>
    <td align="center"><a href="http://about.me/david.nickerson"><img src="https://avatars.githubusercontent.com/u/811244?v=4?s=100" width="100px;" alt=""/><br /><sub><b>David Nickerson</b></sub></a><br /><a href="#ideas-nickerso" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/DeepaMahm"><img src="https://avatars.githubusercontent.com/u/29662579?v=4?s=100" width="100px;" alt=""/><br /><sub><b>DeepaMahm</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/issues?q=author%3ADeepaMahm" title="Bug reports">🐛</a></td>
    <td align="center"><a href="https://github.com/jdieg0"><img src="https://avatars.githubusercontent.com/u/6570972?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Diego</b></sub></a><br /><a href="#tool-jdieg0" title="Tools">🔧</a></td>
    <td align="center"><a href="https://dilawars.me/"><img src="https://avatars.githubusercontent.com/u/895681?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Dilawar Singh</b></sub></a><br /><a href="#tool-dilawar" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/edkerk"><img src="https://avatars.githubusercontent.com/u/7326655?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Eduard Kerkhoven</b></sub></a><br /><a href="#tool-edkerk" title="Tools">🔧</a></td>
    <td align="center"><a href="https://eagmon.github.io/"><img src="https://avatars.githubusercontent.com/u/6809431?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Eran Agmon</b></sub></a><br /><a href="#ideas-eagmon" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/Ermentrout"><img src="https://avatars.githubusercontent.com/u/7952422?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ermentrout</b></sub></a><br /><a href="#tool-Ermentrout" title="Tools">🔧</a></td>
    <td align="center"><a href="https://escher.github.io/"><img src="https://avatars.githubusercontent.com/u/9327950?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Escher</b></sub></a><br /><a href="#tool-escher" title="Tools">🔧</a></td>
    <td align="center"><a href="https://scholar.harvard.edu/fabianfroehlich/home"><img src="https://avatars.githubusercontent.com/u/14923969?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Fabian Fröhlich</b></sub></a><br /><a href="#tool-FFroehlich" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/zhangfengkai"><img src="https://avatars.githubusercontent.com/u/38113699?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Fengkai Zhang</b></sub></a><br /><a href="#tool-zhangfengkai" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/fbergmann"><img src="https://avatars.githubusercontent.com/u/949059?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Frank Bergmann</b></sub></a><br /><a href="#ideas-fbergmann" title="Ideas, Planning, & Feedback">🤔</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/GINsim"><img src="https://avatars.githubusercontent.com/u/32065286?v=4?s=100" width="100px;" alt=""/><br /><sub><b>GINsim</b></sub></a><br /><a href="#tool-GINsim" title="Tools">🔧</a></td>
    <td align="center"><a href="http://gmarupilla.com"><img src="https://avatars.githubusercontent.com/u/53095348?v=4?s=100" width="100px;" alt=""/><br /><sub><b>GMarupilla</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=gmarupilla" title="Code">💻</a></td>
    <td align="center"><a href="http://helikarlab.org/"><img src="https://avatars.githubusercontent.com/u/17307008?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Helikar Lab Personal</b></sub></a><br /><a href="#tool-HelikarLabPersonal" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.sys-bio.org/"><img src="https://avatars.githubusercontent.com/u/1054990?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Herbert Sauro</b></sub></a><br /><a href="#ideas-hsauro" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/hsorby"><img src="https://avatars.githubusercontent.com/u/778048?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Hugh Sorby</b></sub></a><br /><a href="#tool-hsorby" title="Tools">🔧</a></td>
    <td align="center"><a href="http://identifiers.org/"><img src="https://avatars.githubusercontent.com/u/18701545?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Idenfitiers.org</b></sub></a><br /><a href="#data-identifiers-org" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/JanHasenauer"><img src="https://avatars.githubusercontent.com/u/12297214?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jan Hasenauer</b></sub></a><br /><a href="#tool-JanHasenauer" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://bionetgen.org/"><img src="https://avatars.githubusercontent.com/u/8277248?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jim Faeder</b></sub></a><br /><a href="#tool-jrfaeder" title="Tools">🔧</a> <a href="#data-jrfaeder" title="Data">🔣</a></td>
    <td align="center"><a href="http://vcell.org"><img src="https://avatars.githubusercontent.com/u/20616724?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jim Schaff</b></sub></a><br /><a href="#ideas-jcschaff" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/jmrohwer"><img src="https://avatars.githubusercontent.com/u/502289?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Johann Rohwer</b></sub></a><br /><a href="#tool-jmrohwer" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/jhgennari"><img src="https://avatars.githubusercontent.com/u/2684850?v=4?s=100" width="100px;" alt=""/><br /><sub><b>John Gennari</b></sub></a><br /><a href="#ideas-jhgennari" title="Ideas, Planning, & Feedback">🤔</a> <a href="#tool-jhgennari" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/jreadey"><img src="https://avatars.githubusercontent.com/u/7785492?v=4?s=100" width="100px;" alt=""/><br /><sub><b>John Readey</b></sub></a><br /><a href="#tool-jreadey" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/johnsekar"><img src="https://avatars.githubusercontent.com/u/1610689?v=4?s=100" width="100px;" alt=""/><br /><sub><b>John Sekar</b></sub></a><br /><a href="#ideas-johnsekar" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/joncison"><img src="https://avatars.githubusercontent.com/u/1506863?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jon Ison</b></sub></a><br /><a href="#data-joncison" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://www.karrlab.org"><img src="https://avatars.githubusercontent.com/u/2848297?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jonathan Karr</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=jonrkarr" title="Code">💻</a> <a href="https://github.com/biosimulations/biosimulations/commits?author=jonrkarr" title="Documentation">📖</a> <a href="#design-jonrkarr" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/jtcooper10"><img src="https://avatars.githubusercontent.com/u/42880781?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Joshua Cooper</b></sub></a><br /><a href="#tool-jtcooper10" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/starboerg"><img src="https://avatars.githubusercontent.com/u/5522086?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jörn Starruß</b></sub></a><br /><a href="#tool-starboerg" title="Tools">🔧</a></td>
    <td align="center"><a href="http://juergen.pahle.de/"><img src="https://avatars.githubusercontent.com/u/5473011?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jürgen Pahle</b></sub></a><br /><a href="#tool-jpahle" title="Tools">🔧</a></td>
    <td align="center"><a href="https://www.karrlab.org/"><img src="https://avatars.githubusercontent.com/u/13785824?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Karr whole-cell modeling lab</b></sub></a><br /><a href="#ideas-KarrLab" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/0u812"><img src="https://avatars.githubusercontent.com/u/7402146?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Kyle Medley</b></sub></a><br /><a href="#tool-0u812" title="Tools">🔧</a> <a href="#ideas-0u812" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="http://lems.github.io/LEMS"><img src="https://avatars.githubusercontent.com/u/3033237?v=4?s=100" width="100px;" alt=""/><br /><sub><b>LEMS</b></sub></a><br /><a href="#tool-LEMS" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://loicpauleve.name/"><img src="https://avatars.githubusercontent.com/u/228657?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Loïc Paulevé</b></sub></a><br /><a href="#data-pauleve" title="Data">🔣</a> <a href="#tool-pauleve" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/luciansmith"><img src="https://avatars.githubusercontent.com/u/1736150?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Lucian Smith</b></sub></a><br /><a href="#ideas-luciansmith" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/lutzbrusch"><img src="https://avatars.githubusercontent.com/u/13622401?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Lutz Brusch</b></sub></a><br /><a href="#tool-lutzbrusch" title="Tools">🔧</a></td>
    <td align="center"><a href="https://uk.linkedin.com/in/manuelbernal"><img src="https://avatars.githubusercontent.com/u/8855107?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Manuel Bernal Llinares</b></sub></a><br /><a href="#data-mbdebian" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/MarcDinh"><img src="https://avatars.githubusercontent.com/u/50445930?v=4?s=100" width="100px;" alt=""/><br /><sub><b>MarcDinh</b></sub></a><br /><a href="#tool-MarcDinh" title="Tools">🔧</a></td>
    <td align="center"><a href="https://livermetabolism.com/"><img src="https://avatars.githubusercontent.com/u/900538?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Matthias König</b></sub></a><br /><a href="#ideas-matthiaskoenig" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://orcid.org/0000-0002-1509-4981"><img src="https://avatars.githubusercontent.com/u/992660?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Matúš Kalaš</b></sub></a><br /><a href="#data-matuskalas" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/vcellmike"><img src="https://avatars.githubusercontent.com/u/29076280?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michael Blinov</b></sub></a><br /><a href="#ideas-vcellmike" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="http://dumontierlab.com/"><img src="https://avatars.githubusercontent.com/u/993852?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michel Dumontier</b></sub></a><br /><a href="#data-micheldumontier" title="Data">🔣</a></td>
    <td align="center"><a href="http://www.cds.caltech.edu/~mhucka"><img src="https://avatars.githubusercontent.com/u/1450019?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mike Hucka</b></sub></a><br /><a href="#tool-mhucka" title="Tools">🔧</a></td>
    <td align="center"><a href="https://hpc.uchc.edu"><img src="https://avatars.githubusercontent.com/u/400595?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mike Wilson</b></sub></a><br /><a href="#infra-mpw6" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
    <td align="center"><a href="http://modeldb.yale.edu/"><img src="https://avatars.githubusercontent.com/u/38667483?v=4?s=100" width="100px;" alt=""/><br /><sub><b>ModelDB</b></sub></a><br /><a href="#data-ModelDBRepository" title="Data">🔣</a></td>
    <td align="center"><a href="https://unseenbio.com/"><img src="https://avatars.githubusercontent.com/u/135653?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Moritz E. Beber</b></sub></a><br /><a href="#tool-Midnighter" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.nibib.nih.gov/"><img src="https://avatars.githubusercontent.com/u/12418167?v=4?s=100" width="100px;" alt=""/><br /><sub><b>National Institute of Biomedical Imaging and Bioengineering</b></sub></a><br /><a href="#financial-NIBIB" title="Financial">💵</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://nih.gov/"><img src="https://avatars.githubusercontent.com/u/52710462?v=4?s=100" width="100px;" alt=""/><br /><sub><b>National Institutes of Health</b></sub></a><br /><a href="#financial-NIHGOV" title="Financial">💵</a></td>
    <td align="center"><a href="https://nsf.gov/"><img src="https://avatars.githubusercontent.com/u/23663503?v=4?s=100" width="100px;" alt=""/><br /><sub><b>National Science Foundation</b></sub></a><br /><a href="#financial-NSF-open" title="Financial">💵</a></td>
    <td align="center"><a href="https://docs.neuroml.org/"><img src="https://avatars.githubusercontent.com/u/2727519?v=4?s=100" width="100px;" alt=""/><br /><sub><b>NeuroML</b></sub></a><br /><a href="#tool-NeuroML" title="Tools">🔧</a></td>
    <td align="center"><a href="http://neurosimlab.org/"><img src="https://avatars.githubusercontent.com/u/14202113?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Neurosim lab</b></sub></a><br /><a href="#tool-Neurosim-lab" title="Tools">🔧</a></td>
    <td align="center"><a href="https://opencor.ws/"><img src="https://avatars.githubusercontent.com/u/754570?v=4?s=100" width="100px;" alt=""/><br /><sub><b>OpenCOR</b></sub></a><br /><a href="#tool-opencor" title="Tools">🔧</a></td>
    <td align="center"><a href="https://models.physiomeproject.org/"><img src="https://avatars.githubusercontent.com/u/1114929?v=4?s=100" width="100px;" alt=""/><br /><sub><b>PMR2 - the software behind the Auckland Physiome Repository</b></sub></a><br /><a href="#data-PMR2" title="Data">🔣</a></td>
    <td align="center"><a href="http://www.opensourcebrain.org/"><img src="https://avatars.githubusercontent.com/u/1556687?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Padraig Gleeson</b></sub></a><br /><a href="#data-pgleeson" title="Data">🔣</a> <a href="#tool-pgleeson" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/Paytonco"><img src="https://avatars.githubusercontent.com/u/7064808?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Payton Thomas</b></sub></a><br /><a href="#tool-Paytonco" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.comp-sys-bio.org/"><img src="https://avatars.githubusercontent.com/u/2159130?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Pedro Mendes</b></sub></a><br /><a href="#tool-pmendes" title="Tools">🔧</a></td>
    <td align="center"><a href="http://pedromonteiro.org/"><img src="https://avatars.githubusercontent.com/u/2027375?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Pedro T. Monteiro</b></sub></a><br /><a href="#tool-ptgm" title="Tools">🔧</a></td>
    <td align="center"><a href="http://pysces.sourceforge.net/"><img src="https://avatars.githubusercontent.com/u/6103247?v=4?s=100" width="100px;" alt=""/><br /><sub><b>PySCeS: The Python Simulator for Cellular Systems, provides a variety of tools for the analysis of cellular systems</b></sub></a><br /><a href="#tool-PySCeS" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/Ragzz1995"><img src="https://avatars.githubusercontent.com/u/16513966?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Raghul Kannan</b></sub></a><br /><a href="#tool-Ragzz1995" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/rsmsheriff"><img src="https://avatars.githubusercontent.com/u/7849690?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Rahuman Sheriff</b></sub></a><br /><a href="#data-rsmsheriff" title="Data">🔣</a></td>
    <td align="center"><a href="https://raashika03.github.io/rashika.rathi/"><img src="https://avatars.githubusercontent.com/u/45493793?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Rashika Rathi</b></sub></a><br /><a href="#data-raashika03" title="Data">🔣</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://allencell.org/"><img src="https://avatars.githubusercontent.com/u/9079?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ryan Spangler</b></sub></a><br /><a href="#ideas-prismofeverything" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/Ryannjordan"><img src="https://avatars.githubusercontent.com/u/86376602?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ryann Jordan</b></sub></a><br /><a href="#data-Ryannjordan" title="Data">🔣</a></td>
    <td align="center"><a href="http://sbml.org/About"><img src="https://avatars.githubusercontent.com/u/1799692?v=4?s=100" width="100px;" alt=""/><br /><sub><b>SBML Team</b></sub></a><br /><a href="#tool-sbmlteam" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/SED-ML"><img src="https://avatars.githubusercontent.com/u/29736746?v=4?s=100" width="100px;" alt=""/><br /><sub><b>SED-ML</b></sub></a><br /><a href="#tool-SED-ML" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/skeating"><img src="https://avatars.githubusercontent.com/u/1736558?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Sarah Keating</b></sub></a><br /><a href="#tool-skeating" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/shoops"><img src="https://avatars.githubusercontent.com/u/1760522?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Stefan Hoops</b></sub></a><br /><a href="#tool-shoops" title="Tools">🔧</a></td>
    <td align="center"><a href="http://www.smoldyn.org/"><img src="https://avatars.githubusercontent.com/u/33039297?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Steve Andrews</b></sub></a><br /><a href="#data-ssandrews" title="Data">🔣</a> <a href="#tool-ssandrews" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/StochSS"><img src="https://avatars.githubusercontent.com/u/3344600?v=4?s=100" width="100px;" alt=""/><br /><sub><b>StochSS</b></sub></a><br /><a href="#tool-StochSS" title="Tools">🔧</a></td>
    <td align="center"><a href="https://maiage.inrae.fr/en/biosys"><img src="https://avatars.githubusercontent.com/u/32363627?v=4?s=100" width="100px;" alt=""/><br /><sub><b>SysBioINRAe</b></sub></a><br /><a href="#tool-SysBioInra" title="Tools">🔧</a> <a href="#data-SysBioInra" title="Data">🔣</a></td>
    <td align="center"><a href="https://science.vu.nl/en/research/molecular-cell-biology/systems-bioinformatics/index.aspx"><img src="https://avatars.githubusercontent.com/u/12168054?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Systems Biology Lab, Vrije Universiteit Amsterdam</b></sub></a><br /><a href="#tool-SystemsBioinformatics" title="Tools">🔧</a></td>
    <td align="center"><a href="http://systemsbiology.ucsd.edu/"><img src="https://avatars.githubusercontent.com/u/4237829?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Systems Biology Research Group</b></sub></a><br /><a href="#data-SBRG" title="Data">🔣</a></td>
    <td align="center"><a href="https://www.hdfgroup.org/"><img src="https://avatars.githubusercontent.com/u/8572050?v=4?s=100" width="100px;" alt=""/><br /><sub><b>The HDF Group</b></sub></a><br /><a href="#tool-HDFGroup" title="Tools">🔧</a></td>
    <td align="center"><a href="http://neuron.yale.edu/"><img src="https://avatars.githubusercontent.com/u/38567601?v=4?s=100" width="100px;" alt=""/><br /><sub><b>The NEURON Simulator</b></sub></a><br /><a href="#tool-neuronsimulator" title="Tools">🔧</a></td>
    <td align="center"><a href="https://cellml.org/"><img src="https://avatars.githubusercontent.com/u/2141414?v=4?s=100" width="100px;" alt=""/><br /><sub><b>The home of CellML on Github</b></sub></a><br /><a href="#tool-cellml" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/metatoaster"><img src="https://avatars.githubusercontent.com/u/372914?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Tommy Yu</b></sub></a><br /><a href="#data-metatoaster" title="Data">🔣</a></td>
    <td align="center"><a href="https://www.itersdesktop.com/"><img src="https://avatars.githubusercontent.com/u/663341?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Tung Nguyen</b></sub></a><br /><a href="#data-ntung" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/sys-bio"><img src="https://avatars.githubusercontent.com/u/5590646?v=4?s=100" width="100px;" alt=""/><br /><sub><b>UW Sauro Lab</b></sub></a><br /><a href="#tool-sys-bio" title="Tools">🔧</a></td>
    <td align="center"><a href="https://vega.github.io/"><img src="https://avatars.githubusercontent.com/u/11796929?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Vega</b></sub></a><br /><a href="#tool-vega" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/veitveit"><img src="https://avatars.githubusercontent.com/u/15800709?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Veit Schwämmle</b></sub></a><br /><a href="#data-veitveit" title="Data">🔣</a></td>
    <td align="center"><a href="http://vcell.org/"><img src="https://avatars.githubusercontent.com/u/29076025?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Virtual Cell</b></sub></a><br /><a href="#tool-virtualcell" title="Tools">🔧</a> <a href="#data-virtualcell" title="Data">🔣</a></td>
    <td align="center"><a href="http://genome.jouy.inra.fr/~wliebermeis/index_en.html"><img src="https://avatars.githubusercontent.com/u/3976679?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Wolfram Liebermeister</b></sub></a><br /><a href="#tool-liebermeister" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/YinHoon"><img src="https://avatars.githubusercontent.com/u/11270172?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Yin Hoon Chew</b></sub></a><br /><a href="#ideas-YinHoon" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://www.linkedin.com/in/zakandrewking/"><img src="https://avatars.githubusercontent.com/u/1250400?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Zachary A. King</b></sub></a><br /><a href="#tool-zakandrewking" title="Tools">🔧</a> <a href="#data-zakandrewking" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/abulovic"><img src="https://avatars.githubusercontent.com/u/1510530?v=4?s=100" width="100px;" alt=""/><br /><sub><b>abulovic</b></sub></a><br /><a href="#data-abulovic" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/cjmyers"><img src="https://avatars.githubusercontent.com/u/3507191?v=4?s=100" width="100px;" alt=""/><br /><sub><b>cjmyers</b></sub></a><br /><a href="#ideas-cjmyers" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/dczielinski"><img src="https://avatars.githubusercontent.com/u/4442307?v=4?s=100" width="100px;" alt=""/><br /><sub><b>dczielinski</b></sub></a><br /><a href="#tool-dczielinski" title="Tools">🔧</a></td>
    <td align="center"><a href="https://thesustainablevegan.org/"><img src="https://avatars.githubusercontent.com/u/60083977?v=4?s=100" width="100px;" alt=""/><br /><sub><b>freiburgermsu</b></sub></a><br /><a href="https://github.com/biosimulations/biosimulations/commits?author=freiburgermsu" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/jtyurkovich"><img src="https://avatars.githubusercontent.com/u/5396263?v=4?s=100" width="100px;" alt=""/><br /><sub><b>jtyurkovich</b></sub></a><br /><a href="#tool-jtyurkovich" title="Tools">🔧</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://fun.bio.keio.ac.jp/software/libsbmlsim/"><img src="https://avatars.githubusercontent.com/u/16151392?v=4?s=100" width="100px;" alt=""/><br /><sub><b>libsbmlsim</b></sub></a><br /><a href="#tool-libsbmlsim" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/moraru"><img src="https://avatars.githubusercontent.com/u/7397814?v=4?s=100" width="100px;" alt=""/><br /><sub><b>moraru</b></sub></a><br /><a href="#infra-moraru" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#ideas-moraru" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/obodeit"><img src="https://avatars.githubusercontent.com/u/38722594?v=4?s=100" width="100px;" alt=""/><br /><sub><b>obodeit</b></sub></a><br /><a href="#tool-obodeit" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/opencobra"><img src="https://avatars.githubusercontent.com/u/2708410?v=4?s=100" width="100px;" alt=""/><br /><sub><b>openCOBRA</b></sub></a><br /><a href="#tool-opencobra" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/RuleWorld"><img src="https://avatars.githubusercontent.com/u/11491841?v=4?s=100" width="100px;" alt=""/><br /><sub><b>ruleworld</b></sub></a><br /><a href="#tool-RuleWorld" title="Tools">🔧</a> <a href="#data-RuleWorld" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/yexilein"><img src="https://avatars.githubusercontent.com/u/30040612?v=4?s=100" width="100px;" alt=""/><br /><sub><b>yexilein</b></sub></a><br /><a href="#tool-yexilein" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/z-haiman"><img src="https://avatars.githubusercontent.com/u/29131681?v=4?s=100" width="100px;" alt=""/><br /><sub><b>z-haiman</b></sub></a><br /><a href="#tool-z-haiman" title="Tools">🔧</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

A key to the above emojis is available [here](https://allcontributors.org/docs/en/emoji-key).
# Project governance

The project is governed by a committee of core developers who collectively decide the direction of the project with input from the community. A list of the current team members is available [here](../about/team.md).

## Maintainer responsibilities

Maintainers are people who care about the project and are responsible for helping it grow and improve. Maintainers must contribute to the project, as well as help review contributions from the community. Maintainers must also work collaboratively with each other.

Maintainers have write access to this repository. Maintainers can merge their own contributions or contributions from others. Nevertheless, maintainers are encouraged to seek review from each other, particularly for significant changes.

## Becoming a maintainer

To become a maintainer you need to demonstrate the following:

* Participation in project discussions
* Contribution of significant pull requests
* Ability to write quality code, tests, examples, and/or documentation
* Ability to collaborate with the maintainers
* Understanding of the project's goals, organization, and conventions

Prospective maintainers can request maintainer privileges by sending a message to the current maintainers at [info@biosimulations.org](mailto:info@biosimulations.org).

## Project meetings

Maintainers are expected to participate in the project's meetings, which occur online at 11am EST on Thursdays. Other members of the community are also welcome to attend the project's meetings to provide input and feedback on the project. Please contact the maintainers at [info@biosimulations.org](mailto:info@biosimulations.org) for a link to the project meetings.
Fully supported features supported by all BioSimulators tools: 

- Models (`sedml:model`) 
- Model attribute changes (`sedml:changeAttribute`)
- Steady-state (`sedml:steadyState`) and/or uniform timecourse (`sedml:uniformTimeCourse`) simulations
- Algorithms (`sedml:algorithm`) 
- Algorithm parameters (`sedml:algorithmParameter`)
- Tasks for the execution of individual simulations of individual models (`sedml:task`) 
- Data generators for mathematical expressions (`sedml:dataGenerator/@math`) of individual variables (`sedml:dataGenerator`, `sedml:variable`) 
- Reports (`sedml:report`)
- Plots (`sedml:plot2D`)
 
Partially supported, advanced, features supported by some BioSimulators tools:

- Internal model sources (source is a another model in the same SED-ML document)
- Models sourced by identifiers (e.g., `source="biomodels:BIOMD0000001004"`)
- More complex model changes (`sedml:addXML`, `sedml:removeXML`, `sedml:changeXML`, `sedml:computeChange`)
- Repeated tasks (`sedml:repeatedTask`)
- Multi-dimensional reports (`sedml:report`) and plots (`sedml:plot2d`, `sedml:plot3d`)
- 3D plots (`sedml:plot3d`)
- BioNetGen Language  (BNGL)
- CellML
- GINsim Markup Language (GINML, ZGINML)
- NeuroML/Low Entropy Model Specification Language (LEMS)
- Resource Balance Analysis (RBA) XML format
- Smoldyn simulation configuration format
- Systems Biology Markup Language (SBML)
- SBML Flux Balance Constraints (SBML-fbc) package
- SBML Qualitative Models (SBML-qual) package
- SBML Mass Action Stoichiometric Simulation (MASS) schema
To help investigators use Vega to visualize simulation results, BioSimulators-utils provide a command-line program and Python API for converting several types of structural diagrams of models into Vega data visualizations of simulation results. This includes activity network diagrams created with [GINsim](http://ginsim.org/), reaction flux maps created with [Escher](https://escher.github.io/), and [SBGN](https://sbgn.github.io/) process description maps.

* [Interactive tutorial](https://mybinder.org/v2/gh/biosimulators/Biosimulators_tutorials/HEAD?filepath=tutorials/6.%20Generating%20data%20visualizations%20for%20simulation%20results/Converting%20visualizations%20to%20Vega.ipynb) 
* [API documentation](https://docs.biosimulators.org/Biosimulators_utils/source/biosimulators_utils.viz.vega.html)
* [PyPI](https://pypi.org/project/biosimulators-utils/)
* [Source code](https://github.com/biosimulators/Biosimulators_utils)
Several complete examples (COMBINE/OMEX archives with linked SED-ML and Vega files) are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples#compatibility-of-the-example-archives-with-simulation-tools). These examples illustrate how Vega can be combined with a broad range of modeling formalisms, simulation algorithms, model formats, and simulation tools.

* Activity flux diagram created with GINsim with a menu for selecting a perturbation condition (gene knockout) and a slider for scrolling through simulation time
    * Interactive visualization: [synchronous logical updating with SBML-qual, BoolNet](https://biosimulations.org/projects/Yeast-cell-cycle-Irons-J-Theor-Biol-2009#tab=select-viz)
    * COMBINE archive: [synchronous logical updating with SBML-qual](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/sbml-qual/Irons-J-Theor-Biol-2009-yeast-cell-cycle.omex?raw=true)
    * Individual files: [synchronous logical updating with SBML-qual](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/sbml-qual/Irons-J-Theor-Biol-2009-yeast-cell-cycle)
* Line chart with an interactive legend 
    * Interactive visualization: [LSODA with SBML, PySCeS](https://biosimulations.org/projects/Nicotinic-excitation-Edelstein-Biol-Cybern-1996#tab=select-viz)
    * COMBINE archive: [LSODA with SBML](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/sbml-core/Edelstein-Biol-Cybern-1996-Nicotinic-excitation.omex?raw=true)
    * Individual files: [LSODA with SBML](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/sbml-core/Edelstein-Biol-Cybern-1996-Nicotinic-excitation)
* Figure with multiple panels
    * Interactive visualization: [CVODE method with SBML, tellurium](https://biosimulations.org/projects/Morphogenesis-checkpoint-continuous-Ciliberto-J-Cell-Biol-2003#tab=select-viz)
    * COMBINE archive: [CVODE method with SBML](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/sbml-core/Ciliberto-J-Cell-Biol-2003-morphogenesis-checkpoint-continuous.omex?raw=true)
    * Individual files: [CVODE method with SBML](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/sbml-core/Ciliberto-J-Cell-Biol-2003-morphogenesis-checkpoint-continuous)
* Process description map created with SBGN with a slider for scrolling through simulation time 
    * Interactive visualization: [LSODA with SBML, COPASI](https://biosimulations.org/projects/Repressilator-Elowitz-Nature-2000#tab=select-viz)
    * COMBINE archive: [CVODE with CellML](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/cellml/Elowitz-Nature-2000-Repressilator.omex?raw=true) 
    * COMBINE archive: [CVODE with SBML](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/sbml-core/Elowitz-Nature-2000-Repressilator.omex?raw=true)
    * Individual files: [CVODE with CellML](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/cellml/Elowitz-Nature-2000-Repressilator) 
    * Individual files: [CVODE with SBML](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/sbml-core/Elowitz-Nature-2000-Repressilator)
* Reaction flux map created with Escher with interactive metabolites and reactions
    * Interactive visualization: [FBA with SBML-fbc, COBRApy](https://biosimulations.org/projects/Escherichia-coli-core-metabolism-textbook#tab=select-viz)
    * Interactive visualization: [RBA with RBA XML, RBApy](https://biosimulations.org/projects/Escherichia-coli-resource-allocation-Bulovic-Metab-Eng-2019#tab=select-viz)
    * COMBINE archive: [FBA with SBML-fbc](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/sbml-fbc/
    Escherichia-coli-core-metabolism.omex?raw=true)
    * COMBINE archive: [RBA with RBA XML](https://github.com/biosimulators/Biosimulators_test_suite/blob/deploy/examples/rba/Escherichia-coli-K12-WT.omex?raw=true)
    * Individual files: [FBA with SBML-fbc](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/sbml-fbc/
    Escherichia-coli-core-metabolism)
    * Individual files: [RBA with RBA XML](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/rba/Escherichia-coli-K12-WT)
Below are several helpful tools for creating Vega visualizations:

* [Vega examples](https://vega.github.io/vega/examples/)
* [Vega tutorials](https://vega.github.io/vega/tutorials/)
* [Vega documentation](https://vega.github.io/vega/docs/)    
    * [Signals](https://vega.github.io/vega/docs/signals/)
        * [`value` property](https://vega.github.io/vega/docs/signals/#signal-properties)
        * [`bind` property](https://vega.github.io/vega/docs/signals/#bind)
    * [Data sets](https://vega.github.io/vega/docs/data/)
        * [`values` property](https://vega.github.io/vega/docs/data/#properties)
* Software tools for creating Vega visualizations
    - [Altair](https://altair-viz.github.io/) is a Python library which provides methods for generating Vega visualizations. Altair provides similar capabilities to the Matplotib and Seaborn libraries.
    - [Lyra](http://vega.github.io/lyra/) is a interactive graphical program for designing data visualizations.
    - [Vega Editor](https://vega.github.io/editor) is a text editor for Vega documents that continuously renders Vega documents as they are edited. The editor also provides windows for inspecting the values of Vega signals and data sets. The editor is a helpful tool for developing and troubleshooting data visualizations.
    - [Vega-Embed](https://github.com/vega/vega-embed) is a JavaScript program for rendering Vega visualizations inside web pages.
# Live and recorded tutorials

The BioSimulations and BioSimulators Teams periodically present live tutorials. Below is the schedule of upcoming events and recordings of past events.

### HARMONY 2022
We plan to present a tutorial at the HARMONY meeting in Spring 2022.

### COMBINE 2021
Slides from the tutorial are available [here](https://drive.google.com/file/d/1ZituwKeT2jHt4oJDW1DrlanwdDjnWnX6/view?usp=sharing).

### Network modeling summer school 2021
Video of the tutorial on runBioSimulations is available [here](https://www.youtube.com/watch?v=VZGzWBmagcs).

### HARMONY 2021
The tutorial on runBioSimulations and BioSimulators is available from the links below:
* [Video](https://www.youtube.com/watch?v=h0SCn0ZDqq8)
* [Slides](https://drive.google.com/file/d/1Q0X6GNQlT5PfZOcBNVu1B92SSAxJWrnS/view?usp=sharing)
# Guide to using BioSimulations and BioSimulators

This section contains a brief guide for using BioSimulations and BioSimulators to create, publish, and find simulation projects and simulation tools. Please use the menu to the left to navigate through the guide. Below are links to additional tutorials and documentation.

## Live and recorded tutorials
Information about upcoming live tutorials and recordings of past events is available [here](./live-tutorials.md).

## Deeper interactive tutorials on using the REST APIs and Python libraries
Jupyter notebooks with interactive tutorials are available from [Binder](https://tutorial.biosimulators.org/).

* Programmatically introspecting model/simulation files
* Programmatically creating SED-ML files and COMBINE/archives from archetypical model/simulation files
* Programmatically retrieving information about the capabilities of simulation tools
* Programmatically executing simulations in Python

## Example simulation projects
Several example model files, SED-ML files, Vega files, and COMBINE/OMEX archives are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples). This includes examples for all currently supported model languages.

## Documentation for the REST APIs
Links to detailed examples and documentation for the REST APIs is available [here](./api.md).

## Documentation for the BioSimulators simulation tools and core packages
Links to detailed documentation each simulation tool and the core BioSimulators packages are available [here](./biosimulators-packages.md).
# Validating simulation projects (COMBINE/OMEX archives)

BioSimulations provides a tool for validating the contents of a COMBINE/OMEX archive. This tool is available on [runBioSimulations](https://run.biosimulations.org/utils/validate-project). The user can upload a COMBINE/OMEX archive or provide a URL to a COMBINE/OMEX archive. Optionally, the user can select a subset of the archive to validate (e.g., only SED-ML files). The tool will validate the contents of the archive and report any errors. The interface provides a report of the validation results, including validation of the model file, SED-ML file, metadata, archive manifest, and any images presented in the archive.

![validation-tool](./images/validate.png)

!!!warning
    This application can only validate targets for simulation observables for unchanged models. Targets for modified models cannot be validated statically, independently from executing the associated SED-ML document. Simulation tools validate such targets when they execute their parent SED-ML document. For the same reason, this application cannot validate targets for model changes. Simulation tools validate such model change targets during their execution.

This validation can also be accessed programmatically from the [BioSimulations API](https://combine.api.biosimulations.org).
# Programmatically working with BioSimulations and BioSimulators via the REST APIs

APIs are available for submitting and retrieving simulation projects, simulation runs, and simulation tools. Please see the documentation for the [BioSimulations REST API](https://api.biosimulations.org) and [BioSimulators REST API](https://api.biosimulators.org) for more information. The documentation include numerous examples and program for interactively querying the APIs.
# Documentation for the BioSimulators tools and libraries

## Integrated Docker image

A Docker image with a Python environment with most of the validated simulation tools is available at [https://github.com/orgs/biosimulators/packages/container/package/biosimulators](https://github.com/orgs/biosimulators/packages/container/package/biosimulators). An iPython shell for this environment can be launched by installing Docker and running the following commands:

```bash
docker pull ghcr.io/biosimulators/biosimulators
docker run -it --rm ghcr.io/biosimulators/biosimulators
```

This image includes BioSimulators-utils, as well as standardized Python APIs for the simulation tools validated by BioSimulators. Because this image aims to incorporate as many simulation tools as possible within a single Python environment, this image may sometimes lag behind the latest version of each tool.

The Dockerfile for this image is available [here](https://github.com/biosimulators/Biosimulators/blob/dev/Dockerfile).

Information about using the Python APIs in the image is available below.

## Standardized interfaces to simulation tools
Below are links to detailed documentation for the command-line applications and Python APIs for the standardized simulation tools.

* [AMICI](https://docs.biosimulators.org/Biosimulators_AMICI/)
* [Brian 2](https://docs.biosimulators.org/Biosimulators_pyNeuroML/)
* [BioNetGen](https://docs.biosimulators.org/Biosimulators_BioNetGen/)
* [BoolNet](https://docs.biosimulators.org/Biosimulators_BoolNet/)
* [CBMPy](https://docs.biosimulators.org/Biosimulators_CBMPy/)
* [COBRApy](https://docs.biosimulators.org/Biosimulators_COBRApy/)
* [COPASI](https://docs.biosimulators.org/Biosimulators_COPASI/)
* [GillesPy2](https://docs.biosimulators.org/Biosimulators_GillesPy2/)
* [GINsim](https://docs.biosimulators.org/Biosimulators_GINsim/)
* [LibSBMLSim](https://docs.biosimulators.org/Biosimulators_LibSBMLSim/)
* [MASSpy](https://docs.biosimulators.org/Biosimulators_MASSpy/)
* [NetPyNe](https://docs.biosimulators.org/Biosimulators_pyNeuroML/)
* [NEURON](https://docs.biosimulators.org/Biosimulators_pyNeuroML/)
* [OpenCOR](https://docs.biosimulators.org/Biosimulators_OpenCOR/)
* [pyNeuroML](https://docs.biosimulators.org/Biosimulators_pyNeuroML/)
* [PySCeS](https://docs.biosimulators.org/Biosimulators_PySCeS/)
* [RBApy](https://docs.biosimulators.org/Biosimulators_RBApy/)
* [Smoldyn](https://smoldyn.readthedocs.io/en/latest/python/api.html#sed-ml-combine-biosimulators-api)
* [tellurium](https://docs.biosimulators.org/Biosimulators_tellurium/)
* [VCell](https://github.com/virtualcell/vcell)
* [XPP](https://docs.biosimulators.org/Biosimulators_XPP/)

## Template for standardized interfaces to simulation tools
A template repository for creating a standardized interface to a simulation tool is available [here](https://github.com/biosimulators/Biosimulators_simulator_template).

## Core BioSimulators util package
BioSimulators-utils provides (a) several utility command-line programs and a Python API for creating, validating, and executing simulation projects and (b) a Python API for creating standardized interfaces to simulation tools. A tutorial and documentation for the package is available [here](https://docs.biosimulators.org/Biosimulators_utils/).

## Test suite for test simulation tools
A tutorial and documentation for the BioSimulators test suite for simulation tools is available [here](https://docs.biosimulators.org/Biosimulators_test_suite/).
# Creating Vega visualizations of the results of SED-ML files in COMBINE archives

We recommend [Vega](https://vega.github.io/vega/) for data visualizations of simulation results. Vega is a powerful, declarative grammar for describing interactive, two-dimensional data visualizations.

One key feature of Vega is that it modularly captures the graphical marks which comprise visualizations, and how those graphical marks should be painted with data. This feature makes it easy to produce data visualizations for multiple simulation conditions by reusing the same graphical marks with results from multiple simulations. This feature also makes the provenance of data visualizations transparent. As a result, we believe Vega is ideal for collaboration and publication.

Below is a brief tutorial on creating data visualizations with Vega and linking their inputs (Vega data sets) to the outputs of simulations (reports of SED-ML files), as well as several examples. In addition, [BioSimulators-utils](https://github.com/biosimulators/Biosimulators_utils) provides a command-line program and Python API for converting several types of structural diagrams of models into Vega data visualizations of simulation results. See [below](#tools-for-converting-visualizations-of-models-into-vega-data-visualizations) for more information. Links to additional tutorials and documentation for Vega are available [below](#more-information).

## Background SED-ML, COMBINE archive, and Vega concepts

This tutorial focuses on combining Vega data visualizations with simulation experiments described with SED-ML and the COMBINE/OMEX archive format. This requires familiarity with the SED-ML, COMBINE archive, and Vega concepts outlined below. Links to further information about these concepts is also provided below.

* [Simulation Experiment Description Markup Language (SED-ML)](http://sed-ml.org/)
    * Simulations (`sedml:oneStep`, `sedml:steadyState`, `sedml:uniformTimeCourse`)
    * Reports (`sedml:report`)
* [COMBINE/OMEX archive format](https://combinearchive.org/)
    * Manifests which describe the `content` (files) of archives (XML files located at `manifest.xml` in COMBINE archives)
    * The `location` (path) of each `content` (e.g., `path/to/simulation.sedml`)
    * The `format` (media type) of each `content` (e.g., `http://identifiers.org/combine.specifications/sed-ml` for SED-ML, `http://purl.org/NET/mediatypes/application/vnd.vega.v5+json` for Vega)
* [BioSimulations/BioSimulators format for the results of SED-ML reports](../concepts/conventions/simulation-run-reports.md)
* [BioSimulations/BioSimulators conventions for executing COMBINE/OMEX archives with SED-ML files](../concepts/conventions/simulator-interfaces.md)
* [Vega format for data visualizations](https://vega.github.io/vega/)
    * [Signals](https://vega.github.io/vega/docs/signals/)
        * [`value` property](https://vega.github.io/vega/docs/signals/#signal-properties)
        * [`bind` property](https://vega.github.io/vega/docs/signals/#bind)
    * [Data sets](https://vega.github.io/vega/docs/data/)
        * [`values` property](https://vega.github.io/vega/docs/data/#properties)

## Tutorial

The steps below outline how to create Vega data visualizations for simulation results.

1. Create one or more [SED-ML files](https://sed-ml.org/) with one or more reports of the results of one or more simulations.

2. Create a Vega document that captures the desired diagrammatic structure for visualizing these SED-ML reports. We recommend creating Vega documents through tools such as the [Altair Python library](https://altair-viz.github.io/) or the [Lyra graphical editor](https://vega.github.io/lyra/). The [Vega Editor](https://vega.github.io/editor) is also a helpful tool for iteratively editing and troubleshooting Vega visualizations. Furthermore, the [Vega website](https://vega.github.io/vega/) contains extensive examples, tutorials, and documentation for learning Vega.

3. Annotate the Vega signals whose values should be rendered with the values of attributes of simulations or reports of SED-ML documents (e.g., number of a steps of a uniform time course simulation).

    - To set the `value` attribute of a Vega signal equal to the value of an attribute of a simulation or report of a SED-ML document, add a `sedmlUri` key to the signal with a value equal to a list of the location of the SED-ML document, the `id` of the SED-ML simulation or report, and the name of the attribute of the simulation or report (e.g., `['location/of/simulation.sedml', 'simulationId', 'numberOfSteps']`). To indicate that a signal should be rendered with a list of the values of an attribute of multiple simulations or reports, use `SedDocument:*`, `Simulation:*`, or `Report:*` for the SED-ML document location or simulation/report `id` (e.g., `['SedDocument:*', 'Report:*', 'id']` to render a signal with a list of the ids of the all of the reports of all of the SED-ML files in the parent COMBINE/OMEX archive).
    - Similarly, to set the `bind` attribute of a Vega signal equal to the value of an attribute of a simulation or report of a SED-ML document, add a `sedmlUri` key to the `bind` attribute with a value as described above.

4. Annotate the Vega data sets whose values should be rendered with the results of SED-ML reports by adding `sedmlUri` keys to these Vega data sets. The values of these keys should be set as follows to indicate the simulation results that should be linked to each Vega data set:
    - To render a Vega data set with the results of all reports from all of the SED-ML files in the parent COMBINE/OMEX archive, the value of the `sedmlUri` key should be an empty array (i.e. `[]`).
    - To render a Vega data set with the result of a single report from one SED-ML file in the parent COMBINE/OMEX archive, the value of the `sedmlUri` key should be a list of the location of the SED-ML document and the `id` of the report in the document (e.g., `['location/of/simulation.sedml', 'reportId']`).

5. Package the SED-ML and Vega files into a COMBINE/OMEX archive. Include the Vega files in the manifest of the archive with the format `http://purl.org/NET/mediatypes/application/vnd.vega.v5+json`.

## Examples

### Linking visualizations to simulation settings
Below is an example snippet of a Vega document with a regular Vega signal whose value is hard-coded into the Vega document.

```json
{
  "signals": [
    {
      "name": "Regular Vega signal.",
      "value": 10
    }
  ]
}
```

In comparison, below is an example snippet of a Vega document with a signal whose value is dynamically linked to an attribute of a simulation (e.g., number of steps of a time course of a SED-ML file).

```json
{
  "signals": [
    {
      "name": "Vega signal whose value should be rendered with the value of an attribute of a simulation of a SED-ML document.",
      "sedmlUri": [
        "simulation_1.sedml",
        "simulation_1",
        "numberOfSteps"
      ]
    }
  ]
}
```

When BioSimulations renders this Vega document, BioSimulations will retrieve the value of the `numberOfSteps` attribute of the simulation with `id` `simulation_1` of the SED-ML document at location `simulation_1.sedml` (e.g., `10`) and then set the `value` attribute of the signal equal to this value. After this transformation, the value of this signal will be structured as illustrated below.

```json
{
  "signals": [
    {
      "name": "Vega signal whose value should be rendered with the value of an attribute of a simulation of a SED-ML document.",
      "value": 10
    }
  ]
}
```


### Linking visualizations to simulation results

Below is an example snippet of a Vega document with a regular Vega dataset with whose value is hard-coded into the Vega document. Such data sets are useful for capturing the diagrammatic structures of network diagrams, such as the locations of nodes and edges of graphs.

```json
{
  "data": [
    {
      "name": "Regular Vega data set, such as for data for visual marks (e.g., coordinates of nodes and edges of a graph).",
      "values": [
        {
          "id": "A",
          "x": "10",
          "y": "20"
        },
        {
          "id": "B",
          "x": "10",
          "y": "10"
        }
      ]
    }
  ]
}
```

In comparison, below is an example snippet of a Vega document with a data set whose value is dynamically linked to an output of a simulation (report of a SED-ML document). 

```json
{
  "data": [
    {
      "name": "Vega data set whose value should be rendered with the result of a report of a SED-ML document.",
      "sedmlUri": [
        "simulation.sedml",
        "reportId"
      ]
    }
  ]
}
```

When BioSimulations renders this Vega document, BioSimulations will retrieve the value of the report with `id` `reportId` from the SED-ML document at location `simulation.sedml`. The value of this report will be structured as illustrated below.
```json
[
  {
    "id": "t",
    "label": "t",
    "shape": "1,1,11",
    "type": "float64",
    "name": "",
    "values": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  },
  {
    "id": "A",
    "label": "A",
    "shape": "1,1,11",
    "type": "float64",
    "name": "",
    "values": [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
  },
  {
    "id": "B",
    "label": "B",
    "shape": "1,1,11",
    "type": "float64",
    "name": "",
    "values": [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
  }
]
```

After retrieving the value of this report, BioSimulations will set the `values` attribute of the Vega data set to the value of the report. After this transformation, the Vega document will be structured as illustrated below.
```json
{
  "data": [
    {
      "name": "Vega data set whose value should be rendered with the result of a report of a SED-ML document.",
      "values": [
        {
          "id": "t",
          "label": "t",
          "shape": "1,1,11",
          "type": "float64",
          "name": "",
          "values": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        },
        {
          "id": "A",
          "label": "A",
          "shape": "1,1,11",
          "type": "float64",
          "name": "",
          "values": [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        },
        {
          "id": "B",
          "label": "B",
          "shape": "1,1,11",
          "type": "float64",
          "name": "",
          "values": [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
        }
      ]
    }
  ]
}
```

### Complete examples
--8<--
vega-examples.md
--8<--

## Tools for converting visualizations of models into Vega data visualizations
--8<--
vega-converters.md
--8<--

## More information and resources
--8<--
vega-resources.md
--8<--
# Creating simulation projects (COMBINE/OMEX archives with SED-ML files)

## Overview

Authors can publish simulation projects to BioSimulations in four easy steps. Importantly, this workflow supports a broad range of simulations involving a broad range of modeling frameworks, simulation algorithms, model formats, and simulation tools.

1. Create a simulation project (COMBINE/OMEX archive).
    - Create one or more models (e.g., BNGL, CellML, SBML files).
    - Create one or more simulations (SED-ML files).
    - Create zero or more visualizations of the results of the simulations (SED-ML plots or Vega data visualizations).
    - Capture metadata about the project (e.g., taxa, citations, license in an OMEX Metadata RDF-XML file).
2. Execute the project using runBioSimulations.
3. Use runBioSimulations to verify that the project executed as intended.
4. Use runBioSimulations to publish the project.

## Creating projects

BioSimulations uses COMBINE/OMEX archives to encapsulate the files and metadata associated with a project. The guide below summarizes the steps required to create a project. More information about the necessary files, such as the formats, conventions, standards, and organization is provided in the subsequent sections. 

1. Create a directory for the simulation project.

1. Create one or more models: Create files that describe the models that you would like to simulate and save these files to the directory for the project. Alternatively, download a model file from a repository such as BioModels.
    --8<-- "supportedLanguages.md"

1. Create one or more simulations: Create a Simulation Experiment Description Markup Language (SED-ML) file which describes the simulation(s) that you would like to execute and the simulation results that you would like to record.
 
    Currently, all of the validated simulation tools support a subset of SED-ML L1V3 including:

    --8<--
    supportedSedml.md
    --8<--
    
1. Create an OMEX manifest file which describes your model and simulation files: OMEX manifest is a format for describing a collection of files, including the format (e.g., CellML, SBML) of each file. Several example COMBINE archives are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples#compatibility-of-the-example-archives-with-simulation-tools). More information about the OMEX manifest format is available at the COMBINE website.

1. Package your project into a COMBINE archive: A COMBINE archive is a zip file which contains an OMEX manifest file and all of the other files which comprise a project. Several example COMBINE archives are available [here](). More information about the OMEX and COMBINE formats is available at the [COMBINE website](https://combinearchive.org/). The COMBINE website lists several software tools for creating COMBINE archives, including [COMBINE Archive]() and [COMBINE Archive Web](https://cat.bio.informatik.uni-rostock.de/).

Interactive tutorials for creating projects can be found [here](https://tutorial.biosimulators.org). 

## Conventions for simulation projects

To make it easy for investigators to work with a broad range of model formats, modeling frameworks, simulation types, simulation algorithms, and simulation tools, BioSimulations uses existing community standards for creating and describing projects. BioSimulations uses the Simulation Experiment Description Markup Language (SED-ML) to describe [simulation experiments](#sed-ml) and COMBINE/OMEX archive formats for packaging and distributing simulation [projects](#combineomex-archives). Vega is used to describe charts, plots, and [visualizations](#visualizations) of simulation results.  The OMEX-Metadata format is used to describe the [metadata](#metadata) associated with a simulation project. More information about BioSimulation's use of these standards, and how to create and publish a project can be found below.

### SED-ML 

[SED-ML](https://sed-ml.org/) is used to describe simulation experiments. This includes:

- Which models to simulate
- How to modify models to simulate variants, such as alternative initial conditions
- What type of simulations to execute (e.g., steady-state, time course)
- Which algorithms to use (e.g., CVODE, SSA)
- Which observables to record
- How to reduce the recorded values of the observables
- How to plot the observables
- How to export the observables to reports (e.g., CSV, HDF5)

[runBioSimulations](https://run.biosimulations.org/utils/create-project) provides a simple web form for building COMBINE/OMEX archives with SED-ML files from model files (e.g., CellML, SBML). This tool supports all of the model languages supported by BioSimulations. 

The table below provides links to documentation about how to use SED-ML with specific model languages, and provides the URNs and URIs that should be used with SED-ML and OMEX manifests in COMBINE/OMEX archives.

| Language                                                                                                       | EDAM format    | SED-ML URN                    | COMBINE archive specification URI                       | MIME type                 | Extensions        |
| ---------------------------------------------------------------------------------------------------------------|----------------|-------------------------------|---------------------------------------------------------|---------------------------|-------------------|
| [BNGL](https://docs.biosimulators.org/Biosimulators_BioNetGen/tutorial.html)                                   | `format_3972`  | `urn:sedml:language:bngl`     | `http://purl.org/NET/mediatypes/text/bngl+plain`        | `text/bngl+plain`         | `.bngl`           |
| [CellML](http://sed-ml.org/specifications.html)                                                                | `format_3240`  | `urn:sedml:language:cellml`   | `http://identifiers.org/combine.specifications/cellml`  | `application/cellml+xml`  | `.xml`, `.cellml` |
| [(NeuroML)/LEMS](https://docs.neuroml.org/Userdocs/Paths.html)                                                 | `format_9004`  | `urn:sedml:language:lems`     | `http://purl.org/NET/mediatypes/application/lems+xml`   | `application/lems+xml`    | `.xml`            |
| [RBA](https://docs.biosimulators.org/Biosimulators_RBApy/tutorial.html)                                        | `format_9012`  | `urn:sedml:language:rba`      | `http://purl.org/NET/mediatypes/application/rba+zip`    | `application/rba+zip`     | `.zip`            |
| [SBML](http://sed-ml.org/specifications.html)                                                                  | `format_2585`  | `urn:sedml:language:sbml`     | `http://identifiers.org/combine.specifications/sbml`    | `application/sbml+xml`    | `.xml`, `.sbml`   |
| [Smoldyn](https://github.com/ssandrews/Smoldyn/blob/master/Using-Smoldyn-with-SED-ML-COMBINE-BioSimulators.md) | `format_9001`  | `urn:sedml:language:smoldyn`  | `http://purl.org/NET/mediatypes/text/smoldyn+plain`     | `text/smoldyn+plain`      | `.txt`            |
| [XPP](https://docs.biosimulators.org/Biosimulators_XPP/tutorial.html)                                          | `format_9010`  | `urn:sedml:language:xpp`      | `http://purl.org/NET/mediatypes/text/x-xpp`             | `text/x-xpp`              | `.xpp`            |

More detailed information about creating SED-ML files to work with BioSimulations can be found at [here](../concepts/conventions/simulation-experiments.md).

### Visualizations

BioSimulations recommends [Vega](https://vega.github.io/vega/) for data visualizations of simulation results. Vega is a powerful, declarative grammar for describing interactive, two-dimensional data visualizations.

One key feature of Vega is that it modularly captures the graphical marks which comprise visualizations and how those graphical marks should be painted with data. This feature makes it easy to produce data visualizations for multiple simulation conditions by combining the same graphical marks with results from multiple simulations. This features also makes the provenance of data visualizations transparent. As a result, Vega is ideal for collaboration and publication.

Vega files can be included in the COMBINE/OMEX archive of a simulation project to enable visualization on the [project view page](./viewing-projects.md#visualizations).

More information on how to create Vega files that can incorporate data from BioSimulations projects can be found at [here](./creating-vega-visualizations.md).

### Metadata 
BioSimulations uses the [OMEX Metadata](https://co.mbine.org/standards/omex-metadata) specification to capture metadata about the simulation project. BioSimulations has a list of recommended metadata fields that can be used to capture [metadata](./viewing-projects.md#metadata) about a simulation project. Information about how to use these fields can be found [here](../concepts/conventions/simulation-project-metadata.md).

### COMBINE/OMEX archives

BioSimulations uses COMBINE/OMEX archives to bundle the multiple files typically involved in modeling projects into a single archive. These files include, but are not limited to:

- Models (e.g., in CellML, SBML format)
- Simulation experiments (SED-ML files)
- Visualizations for visualizing simulation results (e.g., Vega files)
- Supplementary files, such as data used to calibrate the model
- Metadata about the simulation project (RDF-XML files that follow the OMEX Metadata guidelines)
## Contributing biosimulation tools to BioSimulators

We welcome contributions of additional simulation tools! We encourage developers to submit tools that support BioSimulators' standard containerized command-line interface. Supporting these conventions makes it easier for other investigators to use simulators. However, BioSimulators also accepts simulation tools that are not available via standardized Docker images.

Please follow these steps to contribute a tool to BioSimulators:

1. Annotate the capabilities of your simulation tool (e.g., supported modeling frameworks, simulation algorithms, model formats) using the BioSimulators format for the capabilities of simulation tools.

1. Optionally, build a standardized command-line interface for your simulator. This interface should support the following conventions:

    The command-line interface should support the arguments outlined in BioSimulators' [specifications for command-line interfaces](../concepts/conventions/simulator-interfaces.md) for simulation tools.
    
    - COMBINE/OMEX archives should be used as the format for inputs to your simulator.
    
    - A standard modeling language such as BNGL, CellML, NeuroML, or SBML should be used to describe models.
    
    - SED-ML and the BioSimulators [SED-ML conventions](../concepts/conventions/simulation-experiments.md) should be used to describe simulation experiments.
    
    - Reports of simulation results should be saved according to the [BioSimulators format for reports](../concepts/conventions/simulation-run-reports.md) of simulation results.

    - The process of executing COMBINE/OMEX archives should be logged using [BioSimulators' format for logs](../concepts/conventions/simulation-run-logs.md) of the execution of COMBINE/OMEX archives.
    
    [BioSimulators utils](https://github.com/biosimulators/Biosimulators_utils) provides tools for implementing these conventions. A detailed template for using BioSimulators-utils to build a command-line interface for a simulator is available [here](https://github.com/biosimulators/Biosimulators_simulator_template).

1. Optionally, containerize the command-line interface for your simulator. Such an image will make it easier for others to use your tool. Containerized simulation tools should follow BioSimulators' [conventions for Docker images](../concepts/conventions/simulator-images.md) for simulation tools.

1. Optionally, publish your image to a public repository such as [Docker Hub](https://hub.docker.com/), [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry), or [Quay](https://quay.io/) by executing `docker push { image-url }`. Docker Hub, GitHub Container Registry, and Quay each provide free accounts for public images.

1. [Submit an issue](https://github.com/biosimulators/Biosimulators/issues/new?assignees=&labels=Validate%2Fsubmit+simulator&template=ValidateOrSubmitASimulator.yml&title=%5BSimulation+capabilities%5D%3A+) to the BioSimulators GitHub repository that briefly describes the URL to the specifications of your tool. This will initiate an automated workflow that will validate your simulation tool and either commit your tool to the BioSimulators registry or report problems with your simulation tool that must be addressed. The first version of each simulation tool submitted to the BioSimulators registry will also be manually reviewed by the BioSimulators Team prior to incorporation into the BioSimulators registry.

1. Optionally, set up your continuous integration workflow to automatically push each release to BioSimulators. Within your continuous integration workflow (e.g., CircleCI, GitHub actions, Jenkins, Travis), use the GitHub REST API to automatically create issues that submit each version of your simulator to BioSimulators.

    This requires a GitHub account and a personal access token with the `public_repo` scope. Instructions for creating an access token are available in the [GitHub documentation](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/creating-a-personal-access-token).

    - Endpoint: `http://api.github.com/repos/biosimulators/Biosimulators/issues`
    - Method: `POST`
    - Authentication: `{ github-username }:{ github-person-access-token }`
    - Headers:
        - Accept: `application/vnd.github.v3+json`
    - Body (JSON-encoded)
        - title: `Submit { simulator-name } { simulator-version }`
        - body:

        ```
        ---
        name: { simulator-name }
        version: "{ simulator-version }"
        specificationsUrl: { specifications-url }
        specificationsPatch:
        version: "{ simulator-version }"
        image:
            url: { docker-image-repo-url }:{ simulator-version }
        validateImage: true

        ---
        ```

    To skip the validation of your Docker image, or to submit a simulation tool for which there is no Docker image that provides a BioSimulators-compliant command-line entrypoint, set `validateImage` to `false`.

    Below is an example cURL command to programmatically submit a simulation tool:

    ```bash
    curl \
    -X POST \
    -u jonrkarr:********* \
    -H "Accept: application/vnd.github.v3+json" \
    https://api.github.com/repos/biosimulators/Biosimulators/issues \
    -d '{"labels": ["Validate/submit simulator"], "title": "Submit tellurium 2.1.6", "body": "---\nname: tellurium\nversion: 2.1.6\nspecificationsUrl: https://raw.githubusercontent.com/biosimulators/Biosimulators_tellurium/2.1.6/biosimulators.json\nvalidateImage: true\ncommitSimulator: true\n\n---"}'
    ```
    
    The above is implemented by the Python method `biosimulators_utils.simulator_registry.submit.submit_simulator_to_biosimulators_registry`. See the [documentation](https://docs.biosimulators.org/Biosimulators_utils) for more information.

1. Optionally, also publish the source code for your simulation tool to a repository such as BitBucket, GitHub, or GitLab.
1. Optionally, also publish your simulation tool to a software repository such as CRAN (R), NPM (JavaScript), or PyPi (Python).
1. Optionally, also register your tool with bio.tools. Visit [bio.tools](https://bio.tools/) to submit your tool to their registry of research tools.
1. Optionally, also submit your Dockerfile to BioContainers. BioContainers accepts contributions via pull requests. See the [BioContainers image registry](https://github.com/BioContainers/containers/pulls) for more information.

A sample [continuous integration workflow](https://github.com/biosimulators/Biosimulators_simulator_template/blob/dev/.github/workflows/ci.yml.template) for GitHub Actions is available in the template simulator repository. Instructions for setting up this workflow are in the [README](https://github.com/biosimulators/Biosimulators_simulator_template/blob/dev/README.md).

Each time a commit is pushed to the repository, the workflow executes the following tasks:

1. Clones the repository.

1. Installs the simulator and its dependencies.

1. Lints the code for the simulator.

1. Builds the Docker image for the simulator and tags the image `ghcr.io/<owner>/<repo>/<simulator_id>:<simulator_version` and `ghcr.io/<owner>/<repo>/<simulator_id>:latest`.

1. Runs the unit tests for the simulator and saves the coverage data for the tests.

1. Uploads the coverage results to [Codecov](https://codecov.io/).

1. Compiles the documentation for the simulator.

Each time the repository is tagged (`git tag ...; git push --tags`), the workflow also runs the above tasks. If the above tasks succeed, the workflow executes the following additional tasks:

1. Creates a GitHub release for the tag.
1. Pushes the compiled documentation to the repository (e.g., so it can be served by GitHub pages).
1. Builds the simulator and submits it to a software repository such as PyPI.
1. Pushes the Docker image to the GitHub Container Registry with the above tags.
1. Pushes the simulator to the BioSimulators Registry, via the GitHub API, to create an issue for adding a new version of the simulator to the BioSimulators database. This issue will then automatically use the BioSimulators test suite to validate the simulator and add a new version of the simulator to the database once the simulator passes the test suite.
# Documentation for libraries for publishing primary model repositories to BioSimulations

A portion of BioSimulations is derived from several primary model repositories that focus on models of specific biological systems, specific types of models, specific model formats, and/or specific simulation tools. We aim to automatically publish these repositories to BioSimulations each week. This publication will be implemented by a collection of Python packages and GitHub actions. Links to documentation for these packages will be available below.

| Repository                                                       | Topic                     | Formalism          | Format         | Tools              | Docs for publication to BioSimulations                                         |
|------------------------------------------------------------------|---------------------------|--------------------|----------------|--------------------|---------------------------------------------------------------|
| [BiGG](http://bigg.ucsd.edu/)                                    | Metabolism                | Flux balance       | SBML-fbc       | COBRApy and others | [Docs](https://biosimulations.github.io/biosimulations-bigg/) |
| [BioModels](http://biomodels.net/)                               | Biochemical networks      | Kinetic            | SBML           | Many               |                                                               |
| [GINsim](http://ginsim.org/models_repository)                    | Networks                  | Qualitative        | ZGINML         | GINsim             |                                                               |
| [ModelDB](http://modeldb.science/)                               | Neurophysiology           | Continuous kinetic | ODE and others | XPP and others     |                                                               |
| [Physiome](https://models.physiomeproject.org/)                  | Physiology                | Continuous kinetic | CellML         | OpenCOR and others |                                                               |
| [RBA models](https://github.com/SysBioInra/Bacterial-RBA-models) | Resource allocation       | Resource balance   | RBA XML        | RBAPy              |                                                               |
| [Rule Hub](https://github.com/RuleWorld/RuleHub)                 | Rule-based models         | Discrete kinetic   | BNGL           | BioNetGen          |                                                               |
# Publishing simulation projects (COMBINE/OMEX archives)

## Overview

Authors can publish simulation projects (COMBINE/OMEX archives), their simulation results, and interactive visualizations of their simulation results in three simple steps:

1. Execute the simulation project with runBioSimulations.
2. Review the results of the simulation and their visualizations in runBioSimulations.
3. Use runBioSimulations to publish the run of the project to BioSimulations.

## Executing projects
First, execute your simulation project by uploading it to runBioSimulations and selecting a specific simulation tool. Information about how to use runBioSimulations is available [here](./simulating-projects.md).

## Reviewing the simulation results of projects
Second, once the simulation run has completed, use runBioSimulations to inspect its results. We recommend verifying the following aspects of runs: 

* The COMBINE/OMEX archive includes all relevant files. 
* The simulation behaved as expected and the results are visualized as expected.
* The simulation accurately reproduces biological behavior.
* The simulation project metadata is accurate and complete.
* The thumbnail image displays correctly.

More information about using runBioSimulations to review simulation results is available [here](./viewing-projects.md#visualizations).

## Publishing projects

![share-button](./images/share.png){align=right}

Third, once the simulation run has completed and its results have been reviewed, publish the run to BioSimulations by clicking the "Publish" button on the page for the run. Then select a unique id for the project, consent to BioSimulations' terms, and click the "Publish" button.

## Using BioSimulators to test projects prior to submission to BioSimulations

BioSimulations uses BioSimulators Docker containers to execute simulation projects. This allows authors to run simulation projects locally using the same simulator containers as runBioSimulations. The containerized simulation tools provide a consistent environment and interface that matches that provided on runBioSimulations. This consistency makes it easy for authors to debug problems by enabling authors to use their own machines to interactively run the same simulations as runBioSimulations. For more detailed information and alternative methods to simulating projects see the [simulating projects guide](./simulating-projects.md).

To run a project locally, pull the appropriate BioSimulators Docker image, and then run the simulator as follows:

```bash
docker run ghcr.io/biosimulators/tellurium:2.2.1 -v /path/to/project:/root -v /path/to/output:/root/out -i project.omex -o /root/out
```

!!!note
    Ensure that the directories containing the project.omex file and the output directory are mounted in the container. For more information, see the Docker documentation [here](https://docs.docker.com/storage/bind-mounts/).

## User accounts for owning projects

No login is required to access BioSimulations. However, users must have an account to manage projects. This allows for proper crediting of authors, and it allows authors to manage and edit their projects. These accounts are available for free.

!!!note
    User accounts are under development. If you are interested in an account, please contact us at info@biosimulations.org.

## Privately sharing resources with colleagues, peer reviewers, and editors before publication

![share-button](./images/share.png){align=right}

Before publishing a simulation run, you can share the run privately by providing colleagues the URL of the simulation run on [https://run.biosimulations.org](https://run.biosimulations.org). This link can be retrieved by clicking on the "Share" button in the simulation run view page.
# Frequently asked questions

## Simulation projects (COMBINE/OMEX archive)

**Can I search for projects by wild cards?**

Yes. The [search engine](https::/biosimulations.org/projects) supports wild cards such as `sys*`, `*tems`, and `sy*ems`.

**Can I search for projects by specific attributes?**

Yes. The [search engine](https::/biosimulations.org/projects) supports searching over individual attributes by prepending search queries with the key for the attribute, such as `title:xyz` to search for "xyz" in the title attribute of each project. The key for each attribute is the name of each attribute, lower cased, with spaces replaced by dashes, and without units (e.g., `project-size` for "Project size (MB)"). The table below summarizes the attributes that the search engine currently supports.

| Field                 | Description                                                                                | Units | Key                     | Typical ids                                                                   |
|-----------------------| -------------------------------------------------------------------------------------------|-------|-------------------------|-------------------------------------------------------------------------------|
| Id                    | BioSimulations id                                                                          |       | `id`                    |                                                                               |
| Title                 | Tag line                                                                                   |       | `title`                 |                                                                               |
| Abstract              | Brief summary                                                                              |       | `abstract`              |                                                                               |
| Description           | Extended summary                                                                           |       | `description`           |                                                                               |
| Biology               | Genes, pathways                                                                            |       | `biology`               | [GO](http://geneontology.org/), [UniProt](https://www.uniprot.org/)           |
| Taxa                  |                                                                                            |       | `taxa`                  | [Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)                             |
| Keywords              |                                                                                            |       | `keywords`              |                                                                               |
| Model formats         |                                                                                            |       | `model-formats`         | [EDAM](https://edamontology.org/), [SED-ML URN](https://sed-ml.org/urns.html) |
| Simulation types      |                                                                                            |       | `simulation-types`      |                                                                               |
| Simulation algorithms |                                                                                            |       | `simulation-algorithms` | [KiSAO](https://github.com/SED-ML/KiSAO)                                      |
| Simulation tools      |                                                                                            |       | `simulation-tools`      | [BioSimulators](https://biosimulators.org)                                    |
| Reports               | Report formats                                                                             |       | `reports`               |                                                                               |
| Visualizations        | Visualization formats                                                                      |       | `visualizations`        |                                                                               |
| Project size          | COMBINE archive size                                                                       | MB    | `project-size`          |                                                                               |
| CPUs                  | Requested CPUs                                                                             |       | `cpus`                  |                                                                               |
| Memory                | Requested memory                                                                           | GB    | `memory`                |                                                                               |
| Results size          | Size of outputs                                                                            | MB    | `results-size`          |                                                                               |
| Runtime               |                                                                                            | min   | `runtime`               |                                                                               |
| Simulation provenance | Whether the simulation experiment (SED-ML) was computationally generated from a model file |       | `simulation-provenance` |                                                                               |
| Organizations         |                                                                                            |       | `organizations`         | BioSimulations organization ids                                               |
| Owner                 |                                                                                            |       | `owner`                 | BioSimulations account ids                                                    |
| Creators              | Authors                                                                                    |       | `creators`              | [ORCID](https://orcid.org/)                                                   |
| Contributors          | Curators                                                                                   |       | `contributors`          | [ORCID](https://orcid.org/)                                                   |
| Funders               | Funding agencies                                                                           |       | `funders`               | [Funder Registry](https://www.crossref.org/services/funder-registry/)         |
| Citations             | Publications                                                                               |       | `citations`             | [DOI](https://www.doi.org/), [PubMed](https://pubmed.ncbi.nlm.nih.gov/)       |
| Identifiers           |                                                                                            |       | `identifiers`           | [Identifiers.org](https://identifiers.org/)                                   |
| Additional metadata   |                                                                                            |       | `additional-metadata`   |                                                                               |
| License               |                                                                                            |       | `license`               | [SPDX](https://spdx.org/)                                                     |
| created               | Date archived created                                                                      |       | `created`               |                                                                               |
| published             | Date project published                                                                     |       | `published`             |                                                                               |
| updated               | Date project updated                                                                       |       | `updated`               |                                                                               |

**Which formats for projects does BioSimulations support?**

BioSimulations supports the [COMBINE/OMEX archive](https://combinearchive.org/) format. COMBINE/OMEX archives are zip files that contain an additional manifest file that indicates the format (e.g., CellML, CSV, SBML, SED-ML, PNG, etc.) of each file in the archive. This simple format can capture a broad range of projects. The format also provides authors the flexibility to organize their projects as appropriate.

**Are there any limitations to the COMBINE/OMEX archives and SED-ML files that runBioSimulations can execute?**

runBioSimulations is designed to be able to execute any COMBINE/OMEX archive and the SED-ML files they contain. In practice, runBioSimulations cannot execute every possible COMBINE/OMEX archive and SED-ML file for three main reasons.

First, currently runBioSimulations focuses on the latest version of SED-ML (L1V3) and has limited ability to execute simulation experiments encoded using older versions of SED-ML. Going forward, runBioSimulations will support new versions of SED-ML through new versions of simulation tools submitted to BioSimulators that support these new versions of SED-ML. Because BioSimulators stores old versions of simulation tools, runBioSimulations will also maintain the ability to execute simulations that involve the current version of SED-ML.

Second, runBioSimulations can only execute SED-ML files that involve model formats, modeling frameworks, and simulation algorithms that are supported by at least one of the standardized simulation tools in the BioSimulators registry. Presently, the standardized simulation tools support multiple formats, multiple frameworks, and over 40 algorithms. Going forward, we aim to expand this by working with the community to develop standardized interfaces for additional simulation tools and submitting them to BioSimulators.

Third, runBioSimulations has limited ability to execute SED-ML documents that deviate from the specifications of the COMBINE/OMEX archive and SED-ML formats. In practice, this is the most significant limitation because some simulation tools produce SED-ML files which deviate from the specifications of SED-ML and because most of the existing SED-ML files in model repositories such as BioModels deviate from the specifications of SED-ML (note, we are worked with the BioModels and other teams to correct these issues). Below is a list of known issues which prevent runBioSimulations from executing many SED-ML files.

- Broken model references: The model sources in SED-ML files in some repositories such as BioModels are different from the actual locations of the model files.
- Missing namespace definitions: Most simulation tools do not define URIs for the namespace prefixes used in XPath targets for model variables. Most existing SED-ML files at SED-ML.org and in model repositories also lack definitions of these namespace prefixes.
- Missing attributes: Some simulation tools produce SED-ML files that are missing required attributes.
- Invalid attribute values: Some simulation tools produce SED-ML files that have invalid values of some attributes.
- Non-unique identifiers: Some simulation tools produce SED-ML files that have multiple elements with the same identifier. In such cases, references to these elements cannot be resolved.
- Broken references: Some simulation tools produce SED-ML files that have broken references (e.g., no instance of `sedml:model` has an id attribute equal to the value of the modelReference attribute of a `sedml:task`).
- Invalid XPaths to model variables: Some simulation tools produce SED-ML files that do not have correct XPaths to variables in models. Most frequently, this is because some simulation tools confuse the ids and names of model elements.
- Incorrect KiSAO ids for algorithms: Some simulation tools produce SED-ML files that indicate different algorithms than what the simulation tool actually used to execute the simulation.
- Inconsistent plot axes: Some simulation tools produces SED-ML files where curves in the same plot have mixed (linear and log) x, y, and/or z axes.

runBioSimulations provides two web applications for [validating SED-ML files](https://run.biosimulations.org/utils/validate-simulation) and [COMBINE archives](https://run.biosimulations.org/utils/validate-project). [BioSimulators-utils](https://github.com/biosimulators/biosimulators_utils) provides a command-line program and Python API for validating SED-ML files and COMBINE archives.

We are working with the SED-ML community to clarify the specifications of SED-ML and resolve these inconsistencies. To drive consistency, we also developed the [BioSimulators test suite](https://github.com/biosimulators/Biosimulators_test_suite) for verifying whether simulation tools execute COMBINE/OMEX archives and SED-ML files consistently with the specifications of these formats. In addition, we developed [BioSimulators utils](https://github.com/biosimulators/Biosimulators_utils), a Python package which provides methods for more deeply validating COMBINE/OMEX archives and SED-ML files.

**How can I create a COMBINE/OMEX archive?**

[runBioSimulations](https://run.biosimulations.org) provides a simple web-based tool for creating COMBINE/OMEX archives from model files. [BioSimulators-utils](https://github.com/Biosimulators/Biosimulators_utils) provides a command-line tool and a Python API for creating COMBINE/OMEX archives.

Below are several additional tools for creating SED-ML files and COMBINE/OMEX archives. 

* [CombineArchiveWeb](https://cat.bio.informatik.uni-rostock.de/)
* [COPASI](http://copasi.org/)
* [iBioSim](http://www.async.ece.utah.edu/ibiosim)
* [JWS Online](http://jjj.mib.ac.uk/)
* [tellurium](http://tellurium.analogmachine.org/)
* [VCell](http://vcell.org/)

!!! warning

    Most of these tools are not fully compliant with the SED-ML and COMBINE/OMEX archive standards.


**Can projects include multiple models, simulations, and/or visualizations?**

Yes. Projects can include one or more models, one or more simulations of those models, and zero or more visualizations of the results of those simulations.

**Can I use runBioSimulations to execute simulations involving confidential data (e.g., sensitive patient information)?**

runBioSimulations should not be used for simulations involving confidential data such as information about individual patients. We recommend that investigators who wish to execute simulations involving confidential data use the simulation tools provided by [BioSimulators](https://biosimulators.org/) to execute simulations using their own hardware or contact the [BioSimulators Team](mailto:info@biosimulators.org) to discuss other options.

**Where can I find example COMBINE/OMEX archives?**

BioSimulations provides many archives. In addition, several example archives are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples#compatibility-of-the-example-archives-with-simulation-tools).

**How can I quickly run a sample set of simulations?**

Click [here](https://run.biosimulations.org/simulations?try=1) to load several example simulation projects from the [BioSimulators test suite](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/).

**How can I validate a COMBINE/OMEX archive?**

[runBioSimulations](https://run.biosimulations.org) provides a simple web-based tool for validating COMBINE/OMEX archives. [BioSimulators-utils](https://github.com/Biosimulators/Biosimulators_utils) provides a command-line tool and a Python API for validating COMBINE/OMEX archives.

**How can I find a simulation tool for executing similar simulations on my machine?**

[runBioSimulations](https://run.biosimulations.org) provides a tool for recommending simulation tools for specific simulation projects. In addition, [BioSimulators](https://biosimulators.org) provides tools for searching and filtering for simulation tools that support particular modeling frameworks, model formats, and simulation algorithms.

**How can I execute similar simulations on my machine for further investigation?**

BioSimulations executes projects using simulation tools validated by [BioSimulators](https://biosimulators.org). Each simulation tool is available as a Docker image that provides a consistent command-line interface. In addition, most simulation tools provide consistent Python APIs. More information and tutorials are available from BioSimulators.

**How can I execute a project before publication to BioSimulations?**

BioSimulations executes projects using [runBioSimulations](https://run.biosimulations.org), which uses the simulation tools validated by [BioSimulators](https://biosimulators.org). Investigators can directly use runBioSimulations to execute projects. In addition, each simulation tool is available as a Docker image that provides a consistent command-line interface, and most simulation tools provide consistent Python APIs. These tools can be used to perform the same simulations as BioSimulations when it publishes your project. More information and tutorials are available from BioSimulators and runBioSimulations.

**Where does runBioSimulations execute simulations and store their outputs?**

runBioSimulations executes simulations using a high-performance computing cluster at UConn Health. The runBioSimulations user interface is deployed using Netlify. runBioSimulations' APIs are deployed using Google Cloud. runBioSimulations stores simulation projects (COMBINE/OMEX archives) and their results in MongoDB Atlas.

**Which other tools can I use to execute COMBINE/OMEX archives? Are there any limitations executing archives with other tools?**

runBioSimulations is designed to be compatible with the COMBINE/OMEX archive and SED-ML formats. In principle, this means that other tools such as [COPASI](http://copasi.org/) and [VCell](https://vcell.org/) that support these same formats should be able to execute the same simulations that can be executed with runBioSimulations. In practice, there are two main factors which limit the ability of such tools to execute COMBINE/OMEX archives and SED-ML files.

First, each tool can only execute simulations involving a limited number of model formats, modeling frameworks, and simulation algorithms.

Second, some simulation tools only support a subset of SED-ML. For example, some tools do not support model resolution via URI fragments, some do not support model changes, some do not support repeated tasks, and some do not support 3D plots. Going forward, we hope to use BioSimulators to keep finer-grained track of the capabilities of each tool.

Third, some of the inconsistencies in the implementation of COMBINE and SED-ML described above prevent some simulation tools from consistently executing some COMBINE/OMEX archives and SED-ML files. In particular, some tools cannot execute COMBINE/OMEX archives that organize files into subdirectories.

**Is there a limit to the size of simulation projects that can be published through BioSimulations?**

Simulation projects (COMBINE/OMEX archives) are limited to 1 GB each. Furthermore, uploads are limited to 64 MB; larger archives must be submitted via URLs. In addition, simulation results are currently limited to 5 TB per project. Please contact [runBioSimulations Team](mailto:info@biosimulations.org) to arrange larger projects.

**How much resources are available to each simulation?**

runBioSimulations currently allows users to request up the following amounts of resources for each simulation project:

- Projects: 1 GB
- Cores: 24
- RAM: 192 GB
- Time: 20 days
- Results: 5 TB (The produced HDF5 file and zip archive with HDF5 and plots each must be less than 5 TB)

These limits are easily configurable. Please contact [runBioSimulations Team](mailto:info@biosimulations.org) to arrange additional resources.

**How much resources are available to runBioSimulations?**

runBioSimulations has access to over 50 TFLOPs from 2,168 CPU cores, 11 TB RAM, and 8 PB of storage. Of this, 48 CPU cores and 384 GB RAM are dedicated to runBioSimulations. The remainder is shared with other users of the HPC system at UConn Health. More information about the infrastructure available to runBioSimulations is available from the [UConn Health HPC facility](https://health.uconn.edu/high-performance-computing/resources/).

**How many users can use runBioSimulations simultaneously?**

runBioSimulations is architected to accomodate many simultaneous users. There is no specific limit on the number of users that can use runBioSimulations simultaneously. If users submit many simulations simultaneously, their simulations will be queued and processed in the order in which they were submitted.

The runBioSimulations Team strives to make runBioSimulations' resources are fairly available to all users. If necessary, the runBioSimulations will implement a job priority policy to ensure this.

**Is there a limit to the size of simulation results that runBioSimulations can generate?**

The total size of results output is limited to 5 TB

In addition, the format that runBioSimulations relies on for reports is limited to 32 dimensions. For a non-spatial time course, this means that runBioSimulations can accomodate up to 15 layers of nested repeated tasks. For a 3D spatial time course, this means that runBioSimulations can accomodate up to 13 layers of nested repeated tasks.

The total size of all plots generated from an archive is limited to 5 TB.

**How does runBioSimulations store a list of my simulations?**

The list of simulations that you have submitted is only stored in your local browser. The runBioSimulations server does not maintain user accounts. Unless you optionally provided your email address, the runBioSimulations server has no record of which simulations you requested. As a result, if you clear your browser's cache, you will lose the list of your simulations, and it will not be possible to re-create this list.

**How long does runBioSimulations store simulations?**

runBioSimulations stores simulations and their results permanently. For special cases, contact the [BioSimulators Team](mailto:infor@biosimulators.org) to request deleting simulations and results.

**How can I share projects privately with colleagues and peer reviewers without publishing them?**

[runBioSimulations](https://run.biosimulations.org) provides a unique link for each project. These links can be shared with colleagues, peer reviewers, and editors. These links are not publicly advertised.

**How can I use BioSimulations in conjunction with journal articles?**

We recommend embedding hyperlinks to interactive versions of static figures in figure captions, availability sections, and/or supplementary materials. During peer review, private runBioSimulations hyperlinks can be used as described above. We recommend using Identifiers.org hyperlinks (<code>https://identifiers.org/biosimulations/{project-id}</code>, <code>https://identifiers.org/runbiosimulations/{run-id}</code>).

**Do I need to create an account to publish a project?**

No. Projects can be published anonymously without an account or registration. However, to be able to later edit a project, you must create an account and use that account to publish the project. Once the project is created, only that account will be able to edit the project.

**How can I edit a project that I published?**

The owner of a project can associate the project with new simulation runs. This can be used to correct mistakes and/or provide improved versions. First, use runBioSimulations to create a simulation run. Second, use BioSimulations' [REST API](https://api.biosimulations.org) to modify the project by replacing the old simulation run associated with the project with the new run. The online documentation for the API includes a simple web interface for using the API.

!!! info

    The runBioSimulations website currently only enables investigators to publish simulation runs anonymously. To be able to edit a project, currently, users must initially publish the project using BioSimulations' [REST API](https://api.biosimulations.org). 

Please contact the BioSimulations Team via [email](mailto:info@biosimulations.org) for additional assistance.

**How can I delete a project that I published?**

To ensure projects remain accessible to the community, authors cannot delete projects once they have been published.

Please contact the BioSimulations Team via [email](mailto:info@biosimulations.org) for additional assistance.

**How long does BioSimulations store projects?**

BioSimulations stores projects permanently, including their source files and simulation results.

**How long does runBioSimulations store simulations?**

runBioSimulations stores simulations and their results permanently. For special cases, contact the [BioSimulators Team](mailto:info@biosimulators.org) to request deleting simulations and results.

**How long does runBioSimulations store COMBINE archives created via the online webform?**

Archives created with the webform are temporarily stored for 1 day.

**Why did my simulation fail? How can I get my simulation to succeed?**

There are several reasons why simulations can fail including:

- The simulator that you selected is not capable of executing your archive. Because each simulator supports different modeling frameworks, simulation algorithms, and model formats, any given archive is only compatible with a subset of simulation tools. BioSimulators describes the modeling frameworks, simulation algorithms, and model formats that each simulation tool supports. We recommend using BioSimulators to determine which simulation tools are compatible with your archive. Note, BioSimulators does not capture every detail about the capabilities of each simulation tool, in part, because BioSimulators relies on the community to curate simulation tools. For example, BioSimulators has limited information about which SBML elements (e.g., delays, events) most simulation tools support. As a result, to determine which tools are compatible with your archive, it may also be necessary to read the documentation for several potentially compatible simulation tools.
- Your archive or one of the models or simulations inside your archive is invalid. In particular, because many modeling tools are just beginning to support SED-ML, and some tools do not yet produce valid SED-ML files. We recommend creating simulations with tools that faithfully support SED-ML such as tellurium and VCell.
- Your archive describes one or more simulations that can't be solved. For example, many algorithms may not be able to solve a stiff model up to the desired tolerance in the specified number of execution steps. In this case, we recommend using an alternative algorithm such as CVODE or LSODA for continuous kinetic models or restructuring your model.
- Your archive generates very large outputs. runBioSimulations is architected to support arbitrarily large models and simulations. However, because runBioSimulations hasn't yet been hardened from years of use, users may still discover bugs in runBioSimulations. In this case, please help us improve runBioSimulations by using [GitHub issues](https://github.com/biosimulators/Biosimulators/issues/new/choose) to report problems to the BioSimulators team.

**How can I execute the same simulators as runBioSimulations on my own machine?**
runBioSimulations uses the BioSimulators collection of standardized simulation software tools. BioSimulators containers standardized Docker images for each simulation tool that make it easy to execute simulations. These Docker images are easy to install and run on your own machine. The images can be used on top of any operating system. Please see [https://biosimulators.org](https://biosimulators.org) for more information about how to install and run these images.

Most of the standardized simulation tools in BioSimulators also provide standardized Python APIs. These APIs provide additional flexibility such as combining simulation tools together. A single Docker image with most of the Python APIs is also available. Please see [https://biosimulators.org](https://biosimulators.org) for more information.

**How is each project licensed?**

Each project is provided under the license chosen by its authors. These licenses are displayed on the landing pages for each project. We encourage authors to contribute projects under licenses that permit their reuse, such as the [Creative Commons Universal License (CC0)](https://creativecommons.org/share-your-work/public-domain/cc0/).

**How can I create a badge to embed a project into my website?**

We recommend using Shields.io to generate badges for projects. For example, `https://img.shields.io/badge/BioSimulations-published-green` can be used to generate the badge below.

![Badge](https://img.shields.io/badge/BioSimulations-published-green)

**How can I embed capabilities to create COMBINE/OMEX archives into my website?**

Developers can use runBioSimulations to provide their users capabilities to execute their simulations. Developers can achieve this simply by adding hyperlinks to the create simulation project page, [https://run.biosimulations.org/utils/create-project](https://run.biosimulations.org/utils/create-project).

This page supports several query arguments:

- `modelUrl`: URL for a model file that will be configured in a COMBINE archive. This argument instructs the web form to prefill the model file input with this URL.
- `modelFormat`: EDAM id of the format of the models to execute (e.g., `format_2585` for SBML). This argument instructs the web form select this format.
- `modelingFramework`: SBO id of the modeling framework of the simulations to execute (e.g., `SBO_0000293` for continuous kinetic framework). This argument instructs the web form to select this modeling framework.
- `simulationType`: Name of the type of simulation to create (`OneStep`, `SteadyState`, `UniformTimeCourse`). This argument instructs the web form to select this simulation type.
- `simulationAlgorithm`: KiSAO id of the simulation algorithm to execute (e.g., KISAO_0000019 for CVODE). This argument instructs the web form to select this algorithm.

For example, the URL `https://run.biosimulations.org/run?modelUrl=https%3A%2F%2Fwww.ebi.ac.uk%2Fbiomodels%2Fmodel%2Fdownload%2FBIOMD0000000878.4%3Ffilename%3DLenbury2001.xml&modelFormat=format_2585&modelingFramework=SBO_0000293` can be used to link to the capability to create a COMBINE/OMEX archive for BioModels entry `BIOMD0000000878`.

**How can I embed execution capabilities for my simulations into my website?**

Developers can use runBioSimulations to provide capabilities of execute simulations to their users, by simply adding a hyperlink to [https://run.biosimulations.org/run](https://run.biosimulations.org/run).

The run simulations page supports several query arguments:

- `projectUrl`: URL for a COMBINE/OMEX archive to simulate. This argument instructs the web form to prefill the COMBINE/OMEX archive input with this URL.
- `simulator`: Id of the recommended simulator for executing the COMBINE/OMEX archive (e.g., `tellurium`). This argument instructs the web form to preselect this simulator.
- `simulatorVersion`: Recommended version of the simulator for executing the COMBINE/OMEX archive (e.g., `2.2.1`). This argument instructs the web form to preselect this version.
- `cpus`: Recommended number of CPU cores needed to execute the COMBINE/OMEX archive. This argument instructs the web form to preset the requested number of CPU cores.
- `memory`: Recommended amount of RAM in GB needed to execute the COMBINE/OMEX archive. This argument instructs the web form to preset the requested amount of RAM.
- `maxTime`: Recommended amount of time in minutes needed to execute the COMBINE/OMEX archive. This argument instructs the web form to preset the requested maximum execution time.
- `runName`: Recommended name for the simulation run. This argument instructs the web form to preset the name of the simulation run.

For example, the URL `https://run.biosimulations.org/run?projectUrl=https%3A%2F%2Fwww.ebi.ac.uk%2Fbiomodels%2Fmodel%2Fdownload%2FBIOMD0000000878` can be used to link to the capability to simulate BioModels entry BIOMD0000000878.

**Does runBioSimulations provide an additional service for lower-latent simulation?**

Yes! An endpoint for lower-latent simulation is available [here](https://combine.api.biosimulations.org/). This endpoint is intended for interactive execution of computationally-cheap simulations.

This endpoint runs simulations more quickly by executing simulations in an active Python environment rather than submitting jobs to an HPC queue that use Singularity images to execute simulations.

This endpoint returns simulation results rather than returning an id for later retrieving simulation results. This endpoint can return the outputs of projects in three ways:

- JSON document: includes the results of each SED-ML report and plot
- HDF file: includes the results of each SED-ML report and plot
- Zip file: by default, includes an HDF5 file with the data for each SED-ML report and plot, a PDF file for each plot, a YAML-formatted log of the simulation run, and any additional files produced by the simulation tool

Note, this endpoint has several limitations:

- This endpoint has not been tested as extensively as the main runBioSimulations queue. Please contact the [runBioSimulations Team](mailto:info@biosimulations.org) if you experience issues using this endpoint.
- Currently, only a limited amount of computational resources are dedicated to this endpoint. As a result, this endpoint may be overloaded by frequent simulation.
- This endpoint has access to a limited amount of CPU, memory, and disk.
- Runs are limited to 30 seconds. Longer runs are terminated with timeout errors.
- Only one version of each simulation tool is available. We aim to provide the latest version of each simulation tool possible within a shared Python environment. In cases of conflicting dependencies, this endpoint will be behind the latest version of one or more tools. We aim to avoid this problem by encouraging simulation software developers to update their tools to use recent versions of common dependencies.
- Some validated simulation tools are not available. Because the endpoint is executed within a single Python environment, simulation tools which cannot be installed into a mutual Python environment cannot be executed through this endpoint. Currently, 90% of validated tools are available through the endpoint.
- The endpoint cannot save simulation runs to the runBioSimulations database for later retrieval or sharing.

**Can I use the runBioSimulations API to develop an interactive web application for running simulations?**

Yes! The low-latent simulation endpoint described above is designed to support interactive web applications. Please contact the [runBioSimulations Team](mailto:info@biosimulations.org) to enable CORS access for your application and discuss the computational resources needed for your application.

## Models (e.g., CellML, SBML)

**Which modeling frameworks does BioSimulations support?**

Currently BioSimulations supports constraint-based ([Flux Balance Analysis (FBA)](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000004) and [Resource Balance Analysis (RBA)](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000692)), [continuous kinetic](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000293) (ordinary differential equations (ODE) and differential-algebraic equations (DAE)), [discrete kinetic](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000295) (e.g., Stochastic Simulation Algorithms (SSA)), [logical](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000234), and various [hybrid](https://www.ebi.ac.uk/ols/ontologies/sbo/terms?iri=http%3A%2F%2Fbiomodels.net%2FSBO%2FSBO_0000681) models, including non-spatial, spatial, population-based, and particle-based models. More information about the available simulation methods is available at [BioSimulators](https://biosimulators.org).

**Which model formats does BioSimulations support?**

Currently BioSimulations supports several languages including the [BioNetGen Language (BNGL)](https://bionetgen.org), [CellML](https://cellml.org), the [GINsim](http://ginsim.org/) Markup Language, [NeuroML](https://neuroml.org/)/[Low Entropy Model Specification Langauge (LEMS)](https://lems.github.io/LEMS/), the [RBA XML format](https://sysbioinra.github.io/RBApy/), the [Systems Biology Markup Language (SBML)](https://sbml.org) including the Flux Balance Constraints and Qualitative Models Packages, the [Smoldyn](http://www.smoldyn.org/) simulation configuration format, and the XPP [ODE](http://www.math.pitt.edu/~bard/xpp/help/xppodes.html) format. 

**Which SBML packages does BioSimulations support?**

Currently, BioSimulations supports the core, [Flux Balance Constraints (fbc)](http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/fbc), and [Qualitative Models (qual)](http://sbml.org/Documents/Specifications/SBML_Level_3/Packages/qual) packages.

**How can I contribute an additional modeling framework or model format?**

BioSimulations is extensible to additional modeling frameworks and formats. The community can extend BioSimulations' capabilities by contributing simulation tools to [BioSimulators](https://biosimulators.org). More information, tutorials, and examples are available from BioSimulators.

**How can I request an additional modeling framework or model format?**

Please [submit an issue](https://github.com/biosimulators/Biosimulators/issues/new/choose) to request support for another modeling framework and/or model format.

## Simulation experiments (SED-ML) 

**Which simulation formats does BioSimulations support?**

BioSimulations supports these [conventions](../concepts/conventions/simulation-experiments.md) for the [Simulation Experiment Description Markup Language (SED-ML)](https://sed-ml.org). SED-ML can be used to describe a broad range of simulations involving a broad range of modeling frameworks, model formats, simulation algorithms, and simulation tools.

**Which versions of SED-ML does BioSimulations support?**

BioSimulations all versions of SED-ML. However, BioSimulations does not yet support the new features added to the latest version (L1V4) released this summer. We aim to support these new features soon.

**Which SED-ML features does BioSimulators support?**

All of the simulation tools support the following features:

- Models and model attribute changes: `sedml:model`, `sedml:changeAttribute`.
- At least one of steady-state, one step, and uniform time course simulations: `sedml:steadyState`, `sedml:oneStep`, or `sedml:uniformTimeCourse`.
- Algorithms and their parameters: `sedml:algorithm`, `sedml:algorithmParameter`.
- Tasks for the execution of individual simulations of individual models: `sedml:task`.
- Data generators for individual variables: `sedml:dataGenerator`
- Report outputs: `sedml:report`.

Some of the simulation tools, such as tellurium, support the full SED-ML specification.

**How can I create a SED-ML file?**

[runBioSimulations](https://run.biosimulations.org) provides a simple web-based tool for creating SED-ML documents. [BioSimulators-utils](https://github.com/Biosimulators/Biosimulators_utils) provides a command-line tool and a Python API for creating SED-ML documents.

**How should SED-ML be used with specific model languages?**

Information about how to use SED-ML with specific model languages is available [here](../concepts/conventions/simulation-experiments.md).

**Where can I obtain example COMBINE/OMEX archives and SED-ML files?**

Several examples files are available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples). These examples are verified to work with runBioSimulations. This page lists that simulation tools that are compatible with each example.

**What tools can I use to create SED-ML files?**

runBioSimulations provides an [online tool](https://run.biosimulations.org/utils/create-project) for building SED-ML files. Similar functionality is available as a command-line program and Python API through the [BioSimulators-utils Python package](https://github.com/biosimulators/Biosimulators_utils).

**How can I validate a SED-ML file?**

[runBioSimulations](https://run.biosimulations.org) provides a simple web-based tool for validating SED-ML documents. [BioSimulators-utils](https://github.com/Biosimulators/Biosimulators_utils) provides a command-line tool and a Python API for validating SED-ML documents.

**Which simulation algorithms does BioSimulations support?**

BioSimulations supports all simulation algorithms that are validated by BioSimulators. Currently, this includes over 60 algorithms. See below for contributing an additional simulation algorithm.

**How can I contribute an additional simulation algorithm?**

BioSimulations is extensible to additional simulation algorithms. The community can extend BioSimulations' capabilities by contributing simulation tools to [BioSimulators](https://biosimulators.org). More information, tutorials, and examples are available from BioSimulators.

**How can I request an additional simulation algorithm?**

Please [submit an issue](https://github.com/biosimulators/Biosimulators/issues/new/choose) to request support for another simulation algorithm.

**Can I upload simulation results?**

To ensure the provenance of simulation results, simulation results can only be generated by BioSimulations.

**Can I execute simulations without publishing them?**

Yes. [runBioSimulations](https://run.biosimulations.org/) provides a web application for executing simulations. runBioSimulations does not require an account or registration. Simulations executed with runBioSimulations are private until shared or published.

## Simulation results

**Which formats does BioSimulators support?**

Data for reports and plots are saved in HDF5 format. See the [specifications for data sets](../concepts/conventions/simulation-run-reports.md) for more information.

**How are results encoded into HDF5 files?**

Within HDF5 files, the results of each report (`sedml:report`) and plot (`sedml:plot2D`, `sedml:plot3D`) are saved to paths that are a combination of the relative path of the parent SED-ML document within the parent COMBINE/OMEX archive and the id of the report or plot.

For reports, the rows of these data tables correspond to the data sets (`sedml:dataSet`) specified in the SED-ML definition of the report. The heading of each row is the label of the corresponding data set. For plots, the rows of these data tables correspond to the data generators (`sedml:dataGenerator`) specified in the SED-ML definition of the curves and surfaces of the plot. The heading of each row is the id of the corresponding data generator.

For time course simulations, the columns of these data tables correspond to the predicted time points. Below is a sample data table for a report.

See the [specifications for data tables](../concepts/conventions/simulation-run-reports.md) for more information and examples.

## Data visualizations of simulation results

**Which data visualization formats does BioSimulations support?**

Currently, BioSimulations supports the [Simulation Experiment Description Markup Language (SED-ML)](https://sed-ml.org) and [Vega](https://vega.github.io/vega/) formats. The SED-ML format supports basic line charts. Vega is a powerful format that can be used to describe a broad range of data visualizations, including interactive visualizations that brush complex diagrams with results from multiple individual simulations.

BioSimulations follows these [conventions](../concepts/conventions/simulation-run-visualizations.md) for incorporating Vega visualizations into COMBINE/OMEX archives for simulation projects. These conventions enable authors to describe how the outputs of simulation experiments (simulation results) should be linked to the inputs (data sets) of Vega data visualizations. Importantly, these conventions enable Vega data visualizations to be reused across multiple simulations, such as with different initial conditions, simulation algorithms, simulation tools, models, or modeling frameworks.

**What types of plots does BioSimulations support?**

The SED-ML format supports basic statistical charts. The Vega format supports a broad range of charts, including all canonical statistical charts such as bar and line charts, as well as domain-specific diagrams such as activity flow diagrams, processing description maps, and reaction flux maps.

**What tools are available for creating Vega data visualizations?**

Several tools are available for creating Vega data visualizations, including the [Altair](https://altair-viz.github.io/) Python package, the [Lyra](http://idl.cs.washington.edu/projects/lyra/app/) visual editor, and the Vega [text editor](https://vega.github.io/editor/).

**How can I convert existing diagrams into Vega data visualizations?**

In addition, [BioSimulators-utils](https://github.com/Biosimulators/Biosimulators_utils) provides tools for creating Vega data visualizations from several model visualization formats including [Escher](https://escher.github.io/) metabolic flux maps, [GINsim](http://ginsim.org/) activity flow diagrams, and [Systems Biology Graphical Notation (SGBN)](https://sbgn.github.io/) process description maps.

**How can BioSimulations render additional visualization formats?**

One way to use BioSimulations with additional formats is to convert them to Vega. This can be achieved by writing scripts to convert alternative formats into Vega.

Additional rendering tools could be incorporated into BioSimulations. Please contact the [BioSimulations Team](https://biosimulations.org) to discuss how to integrate additional tools with SED-ML files, COMBINE/OMEX archives, and BioSimulations.

## Metadata about simulation projects

**What metadata is required to publish a project?**

The minimum metadata required for publication is described [here](../concepts/conventions/simulation-project-metadata.md).

**Which formats for metadata does BioSimulations support?**

BioSimulations supports the [RDF-XML format](https://www.w3.org/TR/rdf-syntax-grammar/) and the predicates described [here](../concepts/conventions/simulation-project-metadata.md).

**What tools can be used to create RDF-XML files?**

Several tools are available for creating RDF-XML files:

* [libOmexMeta](https://sys-bio.github.io/libOmexMeta/): creates RDF-XML from metadata in SBML files
* [OpenLink Structured Data Editor](https://osde.openlinksw.com/)
* [RDF Editor](https://dotnetrdf.org/docs/stable/user_guide/Tools-rdfEditor.html)
* [RDF Studio](http://www.linkeddatatools.com/rdf-studio)

## Simulation software tools

**Which simulation tools does BioSimulations support?**

BioSimulations supports all simulation tools validated by [BioSimulators](https://biosimulators.org). Currently, this includes over 20 simulation tools. BioSimulations is extensible to additional simulation tools. The community can extend BioSimulations' capabilities by contributing simulation tools to BioSimulators. More information, tutorials, and examples are available from BioSimulators.

**Do all of the simulators provide consistent Docker images?**

Only the simulators that have a curation status of 5 stars provide Docker images that support the BioSimulators conventions. Going forward, we aim to encourage the community to provide additional standardized Docker images.

**Do all of the simulators provide consistent Python APIs?**

Only the simulators that have a curation status of 5 stars provide Python APIs images that support the BioSimulators conventions. Going forward, we aim to encourage the community to provide additional standardized Docker images.

**How frequently are these tools updated?**

For most tools, BioSimulators immediately incorporates new versions upon their release. These new versions are then immediately available for running and publishing projects.

**Does BioSimulations support old versions of tools?**

Yes. BioSimulations uses BioSimulators to run simulations. BioSimulators stores each version of each tool so that in the future, if needed, old models can be rerun with old versions of tools.

**How can I contribute an additional simulation tool?**

BioSimulations is extensible to additional simulation tools. The community can extend BioSimulations' capabilities by contributing simulation tools to [BioSimulators](https://biosimulators.org). More information, tutorials, and examples are available from BioSimulators.

**Can I submit a simulation tool to the BioSimulators registry for which there is no Docker image, or for which there is a Docker image that is not consistent with the BioSimulators conventions?**
Yes! We encourage the community to submit all simulation tools, even if they do not support the BioSimulators interface and Docker image structure conventions.

Note, BioSimulators can only validate simulation tools that support these conventions. Consequently, when creating GitHub issues to submit simulation tools that do not support the BioSimulators conventions, set `validateImage` to `false` to decline validation of the Docker image (or lack thereof) for your simulation tool. The validation status of your simulation tool (passed or skipped) will be communicated to BioSimulators users.

**Can I submit a simulation tool to the BioSimulators registry for which there is no Python API, or for which there is a Python API image that is not consistent with the BioSimulators conventions?**
Yes! We encourage the community to submit all simulation tools, even if they do not support the BioSimulators conventions.

**How does BioSimulators manage ownership of BioSimulators entries for simulation tools?**
BioSimulators uses GitHub teams to manage ownership of simulation tools.

When a simulator is first submitted to BioSimulators, BioSimulators creates a GitHub Team to own the simulator and Biosimulators adds the submitter as a maintainer of that team. The team will have the name `@biosimulators/{ simulator-id }`. Once the team is created, the submitter will be able to manage the team and add collaborators through [GitHub](https://github.com/orgs/biosimulators/teams/simulator-developers/teams).

Only members of these teams can submit issues to edit the specifications of tools and post new versions of tools. If you are not a member of the team for your simulator, you will need to request access from the maintainers of the GitHub team for your tool.

**How can I request an additional simulation tool?**

Please [submit an issue](https://github.com/biosimulators/Biosimulators/issues/new/choose) to request support for another simulation tool.

**How is each simulation tool licensed?**

The simulation tools available through BioSimulations are provided under the licenses documented for each tool. Please see [BioSimulators](https://biosimulators.org) for more information.

**How can I use a simulator Docker image to execute a simulation?**

Instructions for using simulator Docker images to execute simulations are available [here](../users/simulating-projects.md).

**How can I use a command-line interface image to execute a simulation?**
Instructions for using command-line interfaces to execute simulations are available [here](../users/simulating-projects.md).

**How can I use a Python API to execute a simulation?**
Instructions for using simulator Python APIs to execute simulations are available [here](../users/simulating-projects.md).

In addition, several interactive tutorials are available from [Binder](https://tutorial.biosimulators.org/).

**Is a Docker image available with all of the Python APIs available?**
Yes! `The ghcr.io/biosimulators/biosimulators` image contains most of the Python APIs (90% as of 2021-09-23). More information about the image is available [here](https://github.com/biosimulators/Biosimulators/pkgs/container/biosimulators).

The Dockerfile for this image is available [here](https://github.com/biosimulators/Biosimulators/blob/dev/Dockerfile). To the extent possible, this Docker image uses Pipenv to manage the Python environment inside the image. The `Pipfile` and `Pipfile.lock` files for this environment are available [here](https://github.com/biosimulators/Biosimulators/tree/dev/Dockerfile-assets).

The goal of this image is to provide the latest mutually-compatible versions of the Python APIs for simulation tools. When simulation tools require conflicting versions of dependencies, this image may not have the latest version of each simulation tool. In such cases, individual simulation tools can be updated by running `pipenv install --system --selective-upgrade ...` or `pip install --upgrade ....`. Note, such upgrading may break the functionality of other tools.

In addition, the `ghcr.io/biosimulators/biosimulators-base` image contains the non-Python dependencies for these the Python APIs. More information about the image is available [here](https://github.com/biosimulators/Biosimulators/pkgs/container/biosimulators-base).

**Can the Docker images be run on high-performance computing systems without root access?**
Yes! The Docker images can be run on high-performance computing (HPC) systems where root access is not available by first converting the images to [Singularity](https://sylabs.io/) images. All of the validated images are compatible with Singularity.

**Can the Docker images for simulators be converted to Singularity images for use for high-performance computing systems?**
Yes! All of the validated images are compatible with Singularity. As part of the validation process, we check that each Docker image can be converted into a Singularity image.

**Do any of the simulation tools need commercial licenses?**
No! Currently, no validated simulation tool requires a commercial license. However, some tools can execute simulations more quickly with commercial libraries such as [Gurobi](https://www.gurobi.com/products/gurobi-optimizer/) and [IBM CPLEX](https://www.ibm.com/analytics/cplex-optimizer). Gurobi and IBM both provide free licenses for academic research.

Currently, runBioSimulations can execute COBRApy and RBApy with Gurobi.

**How can I use a commercial license with a simulation tool on my own machine?**
Currently, COBRApy and RBApy can read license keys for Gurobi through environment variables which start with the prefix GRB_. For example, COBRApy and RBApy can use Gurobi (using keys for Gurobi's [Web License Service for Container Environments Only](https://www.gurobi.com/web-license-service/)) inside Docker containers as illustrated below.

```
docker run -it --rm \
  --env \
    GRB_WLSACCESSID=********-****-****-****-************ \
    GRB_WLSSECRET=********-****-****-****-************ \
    GRB_LICENSEID=****** \
  ghcr.io/biosimulators/cobrapy
```

**Can I combine multiple BioSimulators tools together into hybrid or multi-algorithmic simulations?**

Yes! Multiple BioSimulators tools could be used to co-simulate multiple models (potentially involving multiple model languages) using multiple simulation algorithms. We recommend using the Python APIs because they are more flexibile than the command-line interfaces. Specifically, we recommend using the [Vivarium framework](https://vivarium-collective.github.io/) to combine multiple simulation tools and/or models. Vivarium also provides tooling to use BioSimulators.

**Which fields are available for search over simulation tools?**

The [simulation tool search](https://biosimulators.org/simulators) supports search over the following fields:
- `id`
- `name`
- `latest-version`
- `description`
- `frameworks`
- `algorithms`
- `algorithm-parameters`
- `model-formats`
- `simulation-formats`
- `archive-formats`
- `image`
- `cli`
- `python-api`
- `interfaces`
- `operating-systems`
- `programming-languages`
- `curation-status`
- `license`
- `dependencies`
- `authors`
- `citations`
- `identifiers`
- `funding`
- `updated`
- `more-info-url`

**How can I create a badge for a simulation tool to embed into my website?**

We recommend using Shields.io to generate shields for simulators. Several examples are below.

Indicate that a simulator is registered:
Use `https://img.shields.io/badge/BioSimulators-registered-green` to generate the badge below.
![Latest badge](https://img.shields.io/badge/BioSimulators-registered-green)

Indicate that a simulator is valid:
Use `https://img.shields.io/badge/BioSimulators-valid-green` to generate the badge below.
![Latest badge](https://img.shields.io/badge/BioSimulators-valid-green)

Indicate the latest registered version of a simulator:
Use `https://img.shields.io/badge/BioSimulators-{ latest-version }-green` to generate the badge below.
![Latest badge](https://img.shields.io/badge/BioSimulators-2.2.1-green)

## Primary model repositories

**Which model repositories does BioSimulations integrate models from?**

Currently, BioSimulations incorporates models from the [BiGG](http://bigg.ucsd.edu/), [BioModels](http://biomodels.net/), [ModelDB](http://modeldb.science/), and [Physiome](https://models.physiomeproject.org/) model repositories. BioSimulations pulls updates each week.

**How can I contribute models from a model repository?**

The [API](https://api.biosimulations.org) can be used to programmatically contribute projects. Please contact the BioSimulations Team via [email](mailto:info@biosimulations.org) for additional assistance.

## Docker images/containers

**What is a Docker image?**

A Docker image is a template for creating a Docker container. Typically, Docker images are containers for an operating system and additional programs. Docker images enable developers to encapsulate programs and their dependencies into a format that can be shared through online repositories such as Docker Hub or the GitHub Container Registry and run on top of any operating system.

**What is a Docker container?**

A Docker container is a virtual machine. Typically Docker containers are created by instantiating Docker images.

## API

**How can I programmatically retrieve projects and their simulation results?**

BioSimulations' [REST API](https://api.biosimulations.org) can be used to programmatically retrieve projects and their simulation results. Extensive documentation is available online.

The OpenAPI [specification](https://api.biosimulations.org/openapi.json) for the API can be used to create libraries for the API for a broad range of languages. Multiple tools for generating client libraries are available, including [OpenAPI Generator](https://openapi-generator.tech/) and [Swagger Codegen](https://swagger.io/tools/swagger-codegen/).

**How can I programmatically search projects?**

BioSimulations' [REST API](https://api.biosimulations.org) provides an endpoint which returns summaries of each project. This can be used to programmatically search for projects. We plan provide a SPARQL endpoint to support more flexible querying.

## User accounts and organizations

**Is a user account or registration needed to use BioSimulations?**

No account or registration is needed to browse or publish projects to BioSimulations. However, to be able to later edit a project, users must create an account and use that account to publish the project.

**Does runBioSimulations maintain user accounts?**

No. runBioSimulations does not have user accounts. No account or registration is necessary to use runBioSimulations. Optionally, you can provide your email to receive notification of the completion of your simulations.

**How can I create an organization?**

Organizations are groups of users and/or machine accounts. Currently, organizations are used to indicate projects submitted by primary model repositories such as BioModels and consortia. Please contact the [BioSimulations Team](mailto:info@biosimulations.org) to request an organization. In your request, please include the following information:

- Desired id for the organization (lower case letters, numbers, dashes, and underscores)
- Name of the organization
- URL for more information about the organization
- List of need user accounts
  - Desired id (e.g., john-doe)
  - Name (e.g., John Doe)
  - (Optional) URL for more information about the user
- List of needed machine accounts (tokens for programmatically submitting projects)
  - Desired id (e.g., biomodels-bot)
  - Name (e.g., BioModels bot)
  - (Optional) URL for more information about the account

**How does runBioSimulations store a list of my simulations without user accounts?**

The list of your simulations is stored in your local browser. Unless you provided your email address, the runRioSimulations server does not know which simulations you submitted. As a result, if your clear your browser's cache, you will lose the list of your simulations, and it will not be possible to recover this list.
# Executing simulation projects (COMBINE/OMEX archives)

## Using runBioSimulations tool to execute a simulation in the cloud

[runBioSimulations](https://run.biosimulations.org) is a simple application that uses BioSimulators to execute modeling studies. runBioSimulations also provides a REST API for programmatically executing simulations.

### Submit a simulation run

Please follow these steps to execute a simulation project:

1. Open the [project submission form](https://run.biosimulations.org/run).
1. Select a COMBINE/OMEX file to execute.
1. Select a simulation tool and a specific version of that tool.
1. Enter a name for your project. We recommend choosing a descriptive name that will help you recall the purpose of your project. These names will be particularly helpful if you run multiple projects.
1. Optionally, enter your email address to receive notification when your project has completed and is ready for your analysis.
1. Click the "Run" button. After you click the "Run" button, you will receive a URL where you will be able to view the status of your project and retrieve and visualize its results. If you provided an email address, you will be notified by email when your project has completed. This email will contain the same URL.

### View the status of a simulation run

Please follow these steps to view the status of a simulation run:

1. Open the URL provided after you submitted your run.
1. Click the "Overview" tab to view the status of the run.
1. Once the run has completed, click the "Log" tab to view the console log for the execution of the run.

### Retrieve the results of a simulation run

After your run has completed, please follow these steps to retrieve its results:

1. Open the URL provided after you submitted your run.
1. Click the results icon to download the results of the run as a zip archive. This archive will contain a file for each report specified in each SED-ML file in the COMBINE/OMEX archive for your run.

### Visualize the results of a simulation run

After your project has completed, please follow these steps to visualize its results:

1. Open the URL provided after you submitted your project.
1. Click the "Select chart" tab to open a form for choosing which results to visualize.
1. Use the tab to select a pre-defined chart or design a custom chart. More information about creating visualizations is available [here](./creating-projects.md) and [here](./creating-vega-visualizations.md).
1. Click the "View chart" tab to view the selected chart.

### Programmatically executing projects with the REST API

In addition to our web application, a REST API for executing projects is available at [https://api.biosimulations.org/](https://api.biosimulations.org/). This API supports the same simulation tools as the web interface.

### Running example simulation projects

The runBioSimulation app contains a variety of example simulation projects. Click [here](https://run.biosimulations.org/simulations?try=1) to explore runs of these projects. More information about these examples is available [here](https://github.com/biosimulators/Biosimulators_test_suite/tree/deploy/examples/).

## Using containerized simulation tools 

### Execute a simulation locally

The BioSimulators simulation tools can also be used to execute simulations on your own machine. Please follow these steps to use a containerized simulation tool to execute a modeling study on your own machine.

1. Install the Docker container engine: Detailed instructions for all major operating systems are available from the [Docker website](https://docs.docker.com/get-docker/).
1. Download the simulator(s) that you wish to use: From your console, execute `docker pull ghcr.io/biosimulators/{ simulator-id }` for each simulator that you wish to use. This will download the simulators and install them onto your machine.
1. Use the selected simulator(s) to execute simulations and save their results: Execute the following from your console:

```bash
docker run \
  --tty \
  --rm \
  --mount type=bind,source={ path-to-directory-of-COMBINE-archive },target=/tmp/project,readonly \
  --mount type=bind,source={ path-to-save-results },target=/tmp/results \
  ghcr.io/biosimulators/{ simulator-id } \
    --archive /tmp/project/{ name-of-COMBINE-archive } \
    --out-dir /tmp/results
```
Your COMBINE archive should be located at `path-to-directory-of-COMBINE-archive/name-of-COMBINE-archive`.

The results will be saved to `path-to-save-results`. The data for reports and plots will be saved in Hierarchical Data Format 5 (HDF5) format and plots will be saved in Portable Document Format (PDF) and bundled into a single zip archive. See the [specifications for reports](../concepts/conventions/simulation-run-reports.md) for more information about the format of reports.

For reports, the rows of each data table will represent the data sets (`sedml:dataSet`) outlined in the SED-ML definition of the report. The heading of each row will be the label of the corresponding data set. For plots, the rows of each data table will represent the data generators (`sedml:dataGenerator`) outlined in the SED-ML definition of the plot. The heading of each row will be the id of the corresponding data generator.

Report tables of steady-state simulations will have a single column of the steady-state predictions of each data set. Report tables of one step simulations will have two columns that represent the predicted start and end states of each data set. Report tables of time course simulations will have multiple columns that represent the predicted time course of each data set. Report tables of non-spatial simulations will not have additional dimensions. Report tables of spatial simulations will have additional dimensions that represent the spatial axes of the simulation.

### Execute a simulation in an HPC environment

The BioSimulators simulation tools can also be running in high-performance computing (HPC) environments where root access is not available by first converting the Docker images for the tools into Singularity images.

All of the validated images for simulation tools are compatible with Singularity. As part of the validation process, we check that each Docker image can be converted into a [Singularity image](https://sylabs.io/).

The steps below illustrate how Singularity can be used to execute the simulation tools in HPC environments.

1. Install Singularity: Instructions are available at [https://sylabs.io/docs/](https://sylabs.io/docs/).
1. Pull the Docker image by executing `docker pull ghcr.io/biosimulators/{ id }:{ version }`.
1. Convert the Docker image to a Singularity image by executing `singularity pull { /path/to/save/singularity-image.sif } docker://ghcr.io/biosimulators/{ id }:{ version }`.
1. Run the Singularity image by executing `singularity run { /path/to/save/singularity-image.sif } ....`

## Using a command-line interface for a simulation tool to execute a simulation

The command-line interfaces for simulation tools can also be installed and run locally. Note, this typically requires additional effort beyond using the Docker images because it requires installing the dependencies for simulation tools.

Please follow these steps to use a command-line interface for a simulation tool to execute a modeling study.

1. Install Python >= 3.7.
1. Install pip.
1. Install the dependencies for the simulation tool. Links to installation instructions are available from the pages for each simulation tool.
1. Install the command-line application for the simulation tool. From your console, use pip to install the Python package which provides the command-line application. The names of the Python packages which provide the command-line applications are available from the pages for each simulation tool.
1. Use the command-line program to execute a simulation project and save its results: Execute the following from your console:

```bash
biosimulators-{ simulator-id } \
    --archive { /path/to/COMBINE-archive.omex } \
    --out-dir { /path/to/save/outputs }
```

In the above example, the simulation project is located at `/path/to/COMBINE-archive.omex` and the results will be saved to `/path/to/save/outputs`.


## Using a Python API for a simulation tool to execute a simulation

The Python APIs for simulation tools provide additional flexibility beyond their Docker images and command-line interfaces. However, using these APIs typically requires additional effort beyond using the Docker images because it requires installing the dependencies for simulation tools, as well as some knowledge of the data structures used by BioSimulators.

Please follow these steps to use a Python API for a simulation tool to execute a modeling study.

1. Install Python >= 3.7.
1. Install pip.
1. Install the dependencies for the simulation tool. Links to installation instructions are available from the pages for each simulation tool.
1. Install the Python API for the simulation tool. From your console, use pip to install the Python package which provides the Python API. The names of the Python packages which provide the Python APIs are available from the pages for each simulation tool.
1. Open a Python shell.
1. Import the Python API for the simulation tool. Import the Python module which provides the Python API. The names of the Python modules which provide the Python APIs are available from the pages for each simulation tool.
1. Use the Python API to execute a simulation project and save its results: Execute the following from your Python shell:

```python
import { simulator_module }
archive_filename = '{ /path/to/COMBINE-archive.omex }'
output_dirname = '{ /path/to/save/outputs }'
{ simulator_module }.exec_sedml_docs_in_combine_archive(archive_filename, output_dirname)
```

In the above example, the simulation project is located at `/path/to/COMBINE-archive.omex` and the results will be saved to `/path/to/save/outputs`.

The `ghcr.io/biosimulators/biosimulators` Docker image contains most of the available Python APIs inside a single Python environment. An iPython shell to this environment can be launched by executing the following from your console:

```bash
docker pull ghcr.io/biosimulators/biosimulators
docker run -it --rm ghcr.io/biosimulators/biosimulators
```
# Creating standardized interfaces to biosimulation tools
We welcome contributions of additional simulation tools!

All simulation tools submitted for validation by BioSimulators should support at least one modeling language, SED-ML, KiSAO, and the COMBINE/OMEX format.

- Modeling languages. Containers should support at least a subset of at least one modeling language such as BNGL, CellML, Kappa, NeuroML/LEMS, pharmML, SBML, or Smoldyn.

- SED-ML. Currently, containers should support at least the following features of SED-ML L1V3:
    - Model and model attribute changes: `sedml:model`, `sedml:changeAttribute`.
    - At least one of steady-state, one step, or timecourse simulations: `sedml:steadyState`, `sedml:oneStep`, or `sedml:uniformTimeCourse`.
    - Tasks for the execution of individual simulations of individual models: `sedml:task`.
    - Algorithms and algorithm parameters: `sedml:algorithm`, `sedml:algorithmParameter`.
    - Data generators for individual variables: `sedml:dataGenerator`, `sedml:variable`
    - Report outputs: `sedml:report`.

- KiSAO. SED-ML documents interpreted by containers should use KiSAO terms to indicate algorithms and algorithm parameters. As necessary, create [issues](https://github.com/SED-ML/KiSAO/issues/new?assignees=&labels=New+term&template=request-a-term.md&title=) on the KiSAO repository to request terms for additional algorithms and algorithm parameters.

- COMBINE/OMEX. Containers should support the full COMBINE/OMEX specification.

Please follow the steps below to create a containerized simulation tool that adheres to the BioSimulators conventions. Several examples are available from the [BioSimulators GitHub organization](https://github.com/biosimulations/). A template repository with template Python code for a command-line interface and a template for a Dockerfile is available [here](https://github.com/biosimulators/Biosimulators_simulator_template).

1. Optionally, create a Git repository for your command-line interface and Dockerfile.
1. Implement a BioSimulators-compliant command-line interface to your simulation tool. The interface should accept two keyword arguments:

    - `-i, --archive`: A path to a COMBINE archive that contains descriptions of one or more simulation tasks.
    - `-o, --out-dir`: A path to a directory where the outputs of the simulation tasks should be saved. Data for plots and reports should be saved in HDF5 format (see the [specifications for data sets](../concepts/conventions/simulation-run-reports.md) for more information) and plots should be saved in Portable Document Format (PDF) bundled into a single zip archive. Data for reports and plots should be saved to `{ out-dir }/reports.h5` and plots should be saved to `{ out-dir/plots.zip }`. Within the HDF5 file and the zip file, reports and plots should be saved to paths equal to the relative path of their parent SED-ML documents within the parent COMBINE/OMEX archive and the id of the report/plot.

    For reports, the rows of the data tables should correspond to the data sets (`sedml:dataSet`) specified in the SED-ML definition of the report (e.g., time, specific species). The heading of each row should be the label of the corresponding data set.

    For plots, the rows of the data tables should correspond to the data generators (`sedml:dataGenerator`) specified in the SED-ML definition of each curve and surface of the plot (e.g., time, specific species). The heading of each row should be the id of the corresponding data generator.

    Data tables of steady-state simulations should have a single column of the steady-state predictions of each data set. Data tables of one step simulations should have two columns that represent the predicted start and end states of each data set. Data tables of time course simulations should have multiple columns that represent the predicted time course of each data set. Data tables of non-spatial simulations should not have additional dimensions. Data tables of spatial simulations should have additional dimensions that represent the spatial axes of the simulation.

    See the [specifications for interfaces](../concepts/conventions/simulator-interfaces.md) for more information.

    In addition, we recommend providing optional arguments for reporting help and version information about your simulator:

    - `-h, --help`: This argument should instruct the command-line program to print help information about itself.
    - `-v, --version`: This argument should instruct the command-line program to report version information about itself.

    The easiest way to create a BioSimulators-compliant command-line interface is to create a [BioSimulators-compliant Python API](../concepts/conventions/simulator-interfaces.md#conventions-for-python-apis) and then use methods in [BioSimulators-utils](https://github.com/biosimulators/Biosimulators_utils) to build a command-line interface from this API. Implementing a BioSimulators-compliant Python API primarily entails implementing a single method for executing a single simulation of a single model. Additional information about creating BioSimulators-compliant Python APIs, command-line interfaces, and Docker images, including templates, is available [here](https://github.com/biosimulators/Biosimulators_simulator_template).

    Simulation tools can also utilize two environment variables to obtain information about the environment that runBioSimulations uses to execute simulations.

    * `HPC`: runBioSimulations sets this variable to `1` to indicate that the simulation tool is being executed in an HPC environment.
    * `CPUS`: runBioSimulations sets this variable to the number of CPUs allocated to the job in which the simulation tool is being executed.

1. Create a Dockerfile which describes how to build an image for your simulation tool.
    1. Use the `FROM` directive to choose a base operating system such as Ubuntu.
    1. Use the `RUN` directive to describe how to install your tool and any dependencies. Because Docker images are typically run as root, reserve `/root` for the home directory of the user which executes the image. Similarly, reserve `/tmp` for temporary files that must be created during the execution of the image. Install your simulation tool into a different directory than `/root` and `/tmp` such as `/usr/local/bin`.
    1. Ideally, the simulation tools inside images should be installed from internet sources so that the construction of an image is completely specified by its Dockerfile and, therefore, reproducible and portable. Additional files needed during the building of the image, such as licenses to commercial software, can be copied from a local directory such as `assets/`. These files can then be deleted and squashed out of the final image and injected again when the image is executed.
    1. Set the `ENTRYPOINT` directive to the path to your command-line interface.
    1. Set the `CMD` directive to `[]`.
    1. Use the `ENV` directive to declare all [environment variables](../concepts/conventions/simulator-interfaces.md#environment-variables) that your simulation tool supports.
    1. Do not use the `USER` directive to set the user which will execute the image so that the user can be set at execution time.
    1. Use the `LABEL` directive to provide the metadata about your simulation tool described below. This metadata is also necessary to submit your image to [BioContainers](https://biocontainers.pro/), a broad registry of images for biological research.
    
        Open Containers Initiative labels:    

        - `org.opencontainers.image.title`: Human-readable title of the image.
        
        - `org.opencontainers.image.version`: Version of the software in the image.
        
        - `org.opencontainers.image.revision`: Source control revision identifier of the software in the image.

        - `org.opencontainers.image.description`: Human-readable description of the software in the image.

        - `org.opencontainers.image.url`: URL to find more information about the image.

        - `org.opencontainers.image.documentation`: URL to get documentation about the image.

        - `org.opencontainers.image.source`: URL to get the Dockerfile for building the image.

        - `org.opencontainers.image.authors`: Contact details of the people or organization responsible for the image.

        - `org.opencontainers.image.vendor`: Name of the entity, organization or individual which distributes the image.

        - `org.opencontainers.image.licenses`: SPDX expression that describes the license(s) under which the software in the image is distributed.

        - `org.opencontainers.image.created`: Date and time when the image was built (RFC 3339).        
    
        BioContainers labels:

        - `version`: Version of your image (e.g., `1.0.0`)

        - `software`: Simulation program wrapped into your image (e.g., `BioNetGen`).

        - `software.version`: Version of the simulation program wrapped into your image (e.g., `2.5.0`).

        - `about.summary`: Short description of the simulation program (e.g., `Package for rule-based modeling of complex biochemical systems`).

        - `about.home`: URL for the simulation program (e.g., `https://bionetgen.org/`).

        - `about.documentation`: URL for documentation for the simulation program (e.g., `https://bionetgen.org/`).

        - `about.license_file`: URL for the license for the simulation program (e.g., `https://github.com/RuleWorld/bionetgen/blob/master/LICENSE`).

        - `about.license`: SPDX license id for the license for the simulation program (e.g., `SPDX:MIT`). See SPDX  for a list of licenses and their ids.

        - `about.tags`: Comma-separated list of tags which describe the simulation program (e.g., `rule-based modeling,dynamical simulation,systems biology,BNGL,BioSimulators`). Please include the tag `BioSimulators`.
        
        - `extra.identifiers.biotools`: Optionally, the bio.tools identifier for the simulation program (e.g., `bionetgen`). Visit [bio.tools](https://bio.tools) to request an identifier for your simulation program.

        - `maintainer`: Name and email of the person/team who developed the image (e.g., `Jonathan Karr <karr@mssm.edu>`).
        
    Below is an example of metadata for the BioNetGen image.
    ``` dockerfile
    LABEL \
    org.opencontainers.image.title="BioNetGen" \
    org.opencontainers.image.version="2.5.0" \
    org.opencontainers.image.description="Package for rule-based modeling of complex biochemical systems" \
    org.opencontainers.image.url="https://bionetgen.org/" \
    org.opencontainers.image.documentation="https://bionetgen.org/" \
    org.opencontainers.image.source="https://github.com/biosimulators/biosimulators_bionetgen" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="MIT" \
    \
    base_image="ubuntu:18.04"
    version="1.0.0"
    software="BioNetGen"
    software.version="2.5.0"
    about.summary="Package for rule-based modeling of complex biochemical systems"
    about.home="https://bionetgen.org/"
    about.documentation="https://bionetgen.org/"
    about.license_file="https://github.com/RuleWorld/bionetgen/blob/master/LICENSE"
    about.license="SPDX:MIT"
    about.tags="rule-based modeling,dynamical simulation,systems biology,biochemical networks,BNGL,SED-ML,COMBINE,OMEX,BioSimulators"
    extra.identifiers.biotools="bionetgen"
    maintainer="Jonathan Karr <karr@mssm.edu>"
    ```

1. Use Docker and the Dockerfile to build an image for your simulation tool by executing the following from your console:
    ```bash
    docker build \
      --tag { dockerhub-user-id }/{ simulator-id }:{ simulator-version } \
      --tag { dockerhub-user-id }/{ simulator-id }:latest \
      --file { path-to-Dockerfile } \
    ```

1. Create a JSON-encoded file that specifies the capabilities of your simulation tool. This file should adhere to the schema described [here](../concepts/conventions/simulator-capabilities.md).

    Use SBO, KiSAO, and EDAM to describe the modeling frameworks, simulation algorithms, and modeling formats that your simulation tool supports. As necessary, use their issue trackers to request additional [SBO](https://sourceforge.net/p/sbo/term-request/), [KiSAO](https://github.com/SED-ML/KiSAO/issues/new?assignees=&labels=New+term&template=request-a-term.md&title=), and [EDAM](https://github.com/edamontology/edamontology/issues/new?template=new-format.md) terms. As necessary, also use the [SED-ML issue tracker](https://github.com/SED-ML/sed-ml/issues/new) to request URNs for additional modeling languages.

1. Use the BioSimulators test suite to validate your image and its specifications. This will check that the image provides a BioSimulators-compliant Docker structure and command-line interface, and that it provides all of the simulation algorithms described in its specifications.

    ```bash
    pip install biosimulators-test-suite
    biosimulators-test-suite validate { dockerhub-user-id }/{ image-id } { path-to-specifications.json }
    ```

    The command-line program for the test suite provides several helpful options, such as for executing tasks directly through command-line interfaces and for executing individual test cases. More information is available [here](https://github.com/biosimulators/Biosimulators_test_suite).
# Finding, exploring & reusing published projects

## Finding projects

Published projects can be browsed at [https://biosimulations.org/projects](https://biosimulations.org/projects). Each card presents a project, with a thumbnail and title. Mousing over the thumbnail shows additional details about the project. You can customize the attributes that are displayed. 

![browse-projects](./images/browse.png)

### Selecting attributes

![select-attributes](./images/select.png){align=right}
Clicking on the search icon in the top-right corner of the page opens a menu with an attributes sub-menu. From here, you can select the attributes that are displayed. Selecting a field will add the attributes to the details presented in the project card when you mouse over each thumbnail.

 <!-- The new lines must be followed by two spaces-->
&NewLine;  
&NewLine;  
&NewLine;  
&NewLine;  
&NewLine;  
&NewLine;  

### Searching for projects
![search-projects](./images/search.png){align=right}
Clicking on the the search icon at the top right of the page opens a search box. A search term, such as 'metabolism' can be entered in the search box. By default, the search term is searched against each [attribute](#selecting-attributes) of each project. Optionally, you can restrict the search to specific attributes. For example, if you want to search for projects that have the taxa 'Escherichia coli', you can enter 'taxa:Escherichia coli' in the search box. For attributes with spaces in the name, replace these spaces with "-". For example, the term "last-updated:2020" can be used to search for projects that contain the value "2020" in the ;last updated' attribute. A list of the available search fields is available in the [FAQs](faqs.md).

### Filtering projects

![filter-attributes](./images/filter.png){align=right}
The list of displayed projects can be filtered by the values of their attributes. For each available attribute, a menu of values is presented. Selecting a value will filter the list of projects to include only those with the selected values.

 <!-- The new lines must be followed by two spaces-->
&NewLine;  
&NewLine;  
&NewLine;  
&NewLine;  
&NewLine;  
&NewLine;  

## Exploring projects

Clicking on a project card opens a page with the project details. The "Overview" tab provides metadata about the selected project, as well as information about underlying model and simulation run. The "Select chart" tab allows you to configure visualizations of the simulation results that can then be viewed on the "View chart" tab. The "Files" tab provides downloads for the files of the project.

### Metadata

BioSimulations collects metadata to enable searching, browsing and discovering projects. The metadata includes information about authorship, license, funding and other provenance information. It also includes information about the modelled system, such as the modelled organism, and tags that describe the project.

![project-metadata](./images/metadata.png)

### Visualizations

Projects can be visualized using both predefined and custom visualizations. The "Select chart" tab allows you to select from pre-defined visualizations that were included in the project, including both basic charts described with SED-ML and more complex visualizations described with Vega. Additionally, you can create your own custom visualizations by selecting one of the "Design a chart" options including histograms, heatmaps and lineplots. Selecting a plot time will open an additional menu with configuration options to select the datasets to be plotted.

Once you have configured your visualization, you can view it by clicking on the "View chart" button. The "Export to Vega" button will export the visualization to a [Vega](https://vega.github.io/) specification, which enables greater user customization of the visual. More information for using Vega with BioSimulations in available [here](../concepts/conventions/simulation-run-visualizations.md).

More information about creating SED-ML and Vega visualizations is available [here](./creating-projects.md) and [here](./creating-vega-visualizations.md).

### Simulation runs
![Sidebar screenshot showing simulation run details](./images/sidebar-simulation-run.png){ align=left }
More detailed information about the execution of the project and its results can be viewed by following the links to the runBioSimulations page for the project. The "Logs" tab provides detailed output of the simulation execution, including each individual simulation task and the outputs (reports and plots) produced by the simulation. Each task of the simulation is presented as a collapsible section that can be expanded to show the outputs. Both structured log files and raw output files can be downloaded from the links. 

## Reusing projects

### Creating and executing variants of simulations with runBioSimulations

In addition to this full-featured web application, [runBioSimulations](https://run.biosimulations.org) provides a simpler web application and REST API for executing simulations. runBioSimulations simply enables users to execute COMBINE archives using a variety of simulation tools and generate time series plots of their results. runBioSimulations does not require an account.

### Downloading projects and executing them with your own computers

#### Downloading projects

The models, simulations, and visualizations in BioSimulations can be programmatically obtained using our [REST API](https://api.biosimulations.org). Documentation for the API is available at the same URL.

#### Recommended tools for further exploring simulation projects

BioSimulations provides basic capabilities for reproducing and reusing a wide range of biomodeling projects. For further work, we encourage users to use the domain-specific online platforms, desktop programs, and libraries outlined below. Consistent interfaces to the desktop and library tools below are available from [BioSimulators](https://biosimulators.org), including Docker images, command-line interfaces and Python APIs. More information about obtaining and using these tools is available from [BioSimulators](https://biosimulators.org). 

!!! warning

    While the BioSimulators interfaces to these tools support SED-ML and the COMBINE/OMEX archive format, the primary versions of most of the tools below do not support these formats or do not support them consistently with the specifications of the SED-ML format.

    
| Framework          | Language | Online programs                                | Desktop programs                         | Libraries  |
|--------------------|----------|------------------------------------------------|------------------------------------------|------------|
| Continuous kinetic | BNGL     |                                                | [BioNetGen](https://bionetgen.org/)      | [pyBioNetGen](https://pybionetgen.readthedocs.io/)    |
| Continuous kinetic | CellML   |                                                | [OpenCOR](https://opencor.ws/)           | [OpenCOR](https://opencor.ws/)    |
| Continuous kinetic | NeuroML  |                                                | [NetPyNe](http://www.netpyne.org/), [NEURON](https://neuron.yale.edu/neuron/), [pyNeuroML](https://github.com/NeuroML/pyNeuroML)  | [NetPyNe](http://www.netpyne.org/), [NEURON](https://neuron.yale.edu/neuron/), [pyNeuroML](https://github.com/NeuroML/pyNeuroML)    |
| Continuous kinetic | SBML     | [JWS Online](http://jjj.biochem.sun.ac.za/)    | [BioNetGen](https://bionetgen.org/), [COPASI](http://copasi.org/), [tellurium](http://tellurium.analogmachine.org/), [VCell](https://vcell.org/) | [AMICI](https://amici.readthedocs.io/), [GillesPy2](https://stochss.github.io/GillesPy2/), [libRoadRunner](https://libroadrunner.org/), [LibSBMLSim](http://fun.bio.keio.ac.jp/software/libsbmlsim/), [pyBioNetGen](https://pybionetgen.readthedocs.io/), [PySCeS](http://pysces.sourceforge.net/)   |
| Continuous kinetic | XPP ODE  |                                                | [XPP](http://www.math.pitt.edu/~bard/xpp/xpp.html)        |         |
| Discrete kinetic   | BNGL     |                                                | [BioNetGen](https://bionetgen.org/)      | [pyBioNetGen](https://pybionetgen.readthedocs.io/)    |
| Discrete kinetic   | SBML     | [StochSS](https://stochss.org/)                | [BioNetGen](https://bionetgen.org/), [COPASI](http://copasi.org/), [tellurium](http://tellurium.analogmachine.org/), [VCell](https://vcell.org/) | [GillesPy2](https://stochss.github.io/GillesPy2/), [libRoadRunner](https://libroadrunner.org/), [pyBioNetGen](https://pybionetgen.readthedocs.io/)   |
| Flux balance       | SBML     | [Fluxer](https://fluxer.umbc.edu/)             | [CBMPy](http://cbmpy.sourceforge.net/)             | [CBMPy](http://cbmpy.sourceforge.net/), [COBRApy](https://opencobra.github.io/cobrapy/)        |
| Logical            | GINsim   |                                                | [GINsim](http://ginsim.org/)           |    |
| Logical            | SBML     | [Cell Collective](https://cellcollective.org/) | [GINsim](http://ginsim.org/)           | [BoolNet](https://sysbio.uni-ulm.de/?Software:BoolNet#:~:text=BoolNet%20is%20an%20R%20package,available%20from%20BoolNet's%20CRAN%20page.)   |
| MASS               | SBML     |                                                |                                                    | [MASSpy](https://masspy.readthedocs.io/)        |
| Resource balance   | RBA XML  |                                                |                                                    | [RBApy](https://sysbioinra.github.io/RBApy/)        |
| Spatial discrete   | Smoldyn  |                                                | [Smoldyn](https://www.smoldyn.org/)                | [Smoldyn](https://www.smoldyn.org/)        |

<div class="logos">
<div class="logos-row">
    <a
    href="https://cellcollective.org/"
    rel="noopener" target="_blank"
    title="Cell Collective"
    >
    <img
        class="zoom"
        src="/assets/images/about/partners/cell-collective.png"
    />
    </a>

    <a href="https://fluxer.umbc.edu/" rel="noopener" target="_blank" title="Fluxer">
    <img class="zoom" src="/assets/images/about/partners/fluxer.svg" />
    </a>

    <a
    href="https://jjj.biochem.sun.ac.za/"
    rel="noopener" target="_blank"
    title="JWS Online"
    >
    <img class="zoom" src="/assets/images/about/partners/jws.svg" />
    </a>

    <a href="https://stochss.org/" rel="noopener" target="_blank" title="StochSS">
    <img class="zoom" src="/assets/images/about/partners/stochss.svg" />
    </a>

    <a
    href="https://vivarium-collective.github.io"
    rel="noopener" target="_blank"
    title="Vivarium"
    >
    <img class="zoom" src="/assets/images/about/partners/vivarium.svg" />
    </a>
</div>
</div>
# Finding simulation tools

<script type="application/ld+json">{
  "@context": "https://schema.org",
  "@type": "HowTo",
  "name": "How to build and execute simulation projects",
  "abstract": "Guide to building simulation projects with the Simulation Experiment Description Markup Language (SED-ML) and COMBINE/OMEX archive format, finding simulation tools capable of executing specific projects, and using those tools to execute simulations.",
  "keywords": [
    "computational biology",
    "systems biology",
    "mathematical model",
    "numerical simulation",
    "COMBINE",
    "OMEX",
    "Simulation Experiment Description Markup Language",
    "SED-ML",
    "CellML",
    "Systems Biology Markup Language",
    "SBML",
    "Kinetic Simulation Algorithm Ontology",
    "KiSAO",
    "Hierarchical Data Format",
    "HDF5"
  ],
  "tool": [
    {
      "@type": "HowToTool",
      "name": "BioSimulations",
      "description": "Open registry of biological simulation projects.",
      "url": "https://biosimulations.org"
    },
    {
      "@type": "HowToTool",
      "name": "BioSimulators",
      "description": "Open registry of biological simulation software tools.",
      "url": "https://biosimulators.org"
    },
    {
      "@type": "HowToTool",
      "name": "runBioSimulations",
      "description": "Web application for executing biological simulations.",
      "url": "https://run.biosimulations.org"
    }
  ],
  "step": [
    {
      "@type": "HowToStep",
      "name": "Find or build a simulation project.",
      "text": "Obtain a project from a repository such as BioSimulations or use a tool such as runBioSimulations to encode a simulation experiment into the Simulation Experiment Markup Language (SED-ML) and COMBINE/OMEX archive format."
    },
    {
      "@type": "HowToStep",
      "name": "Find a simulation tool that has the capabilities to execute the simulation project.",
      "text": "Use the BioSimulators registry or the runBioSimulations simulator recommendation tool to find a simulation tool that supports the model formats, modeling frameworks, and simulation algorithms required for the project."
    },
    {
      "@type": "HowToStep",
      "name": "Obtain the simulation tool.",
      "text": "Navigate your browser to runBioSimulations, or use Docker or pip to install the simulation tool onto your own machine."
    },
    {
      "@type": "HowToStep",
      "name": "Use the simulation tool to execute the simulation and export its outputs.",
      "text": "Follow the online instructions for runBioSimulations or use the Docker image or Python package to execute the project."
    },
    {
      "@type": "HowToStep",
      "name": "Visualize and analyze the simulation results.",
      "text": "View the generated visualizations in the runBioSimulations web application or the generated PDF files."
    }
  ],
  "educationalLevel": "advanced",
  "estimatedCost": {
    "@type": "MonetaryAmount",
    "value": 0,
    "currency": "USD"
  }
}</script><script type="application/ld+json">{
  "@context": "https://schema.org",
  "@type": "HowTo",
  "name": "How to build and execute simulation projects",
  "abstract": "Guide to building simulation projects with the Simulation Experiment Description Markup Language (SED-ML) and COMBINE/OMEX archive format, finding simulation tools capable of executing specific projects, and using those tools to execute simulations.",
  "keywords": [
    "computational biology",
    "systems biology",
    "mathematical model",
    "numerical simulation",
    "COMBINE",
    "OMEX",
    "Simulation Experiment Description Markup Language",
    "SED-ML",
    "CellML",
    "Systems Biology Markup Language",
    "SBML",
    "Kinetic Simulation Algorithm Ontology",
    "KiSAO",
    "Hierarchical Data Format",
    "HDF5"
  ],
  "tool": [
    {
      "@type": "HowToTool",
      "name": "BioSimulations",
      "description": "Open registry of biological simulation projects.",
      "url": "https://biosimulations.org"
    },
    {
      "@type": "HowToTool",
      "name": "BioSimulators",
      "description": "Open registry of biological simulation software tools.",
      "url": "https://biosimulators.org"
    },
    {
      "@type": "HowToTool",
      "name": "runBioSimulations",
      "description": "Web application for executing biological simulations.",
      "url": "https://run.biosimulations.org"
    }
  ],
  "step": [
    {
      "@type": "HowToStep",
      "name": "Find or build a simulation project.",
      "text": "Obtain a project from a repository such as BioSimulations or use a tool such as runBioSimulations to encode a simulation experiment into the Simulation Experiment Markup Language (SED-ML) and COMBINE/OMEX archive format."
    },
    {
      "@type": "HowToStep",
      "name": "Find a simulation tool that has the capabilities to execute the simulation project.",
      "text": "Use the BioSimulators registry or the runBioSimulations simulator recommendation tool to find a simulation tool that supports the model formats, modeling frameworks, and simulation algorithms required for the project."
    },
    {
      "@type": "HowToStep",
      "name": "Obtain the simulation tool.",
      "text": "Navigate your browser to runBioSimulations, or use Docker or pip to install the simulation tool onto your own machine."
    },
    {
      "@type": "HowToStep",
      "name": "Use the simulation tool to execute the simulation and export its outputs.",
      "text": "Follow the online instructions for runBioSimulations or use the Docker image or Python package to execute the project."
    },
    {
      "@type": "HowToStep",
      "name": "Visualize and analyze the simulation results.",
      "text": "View the generated visualizations in the runBioSimulations web application or the generated PDF files."
    }
  ],
  "educationalLevel": "advanced",
  "estimatedCost": {
    "@type": "MonetaryAmount",
    "value": 0,
    "currency": "USD"
  }
}</script>

## Finding a simulation tool that is capable of executing a modeling project
The simulation tools in the BioSimulators collection support different model formats, modeling frameworks, simulation types, simulation algorithms, and observables. To find a simulation tool, use the following steps:

1. Determine the format and framework of your model, and identify the simulation type and algorithm that you would like to execute. 
1. Browse the [simulators](https://biosimulators.org/simulators) to identify a tool which can execute your project. 

Alternatively, runBioSimulations provides a [utility](https://run.biosimulations.org/utils/suggest-simulator) to recommend a simulation tool for specific combinations of model formats, modeling frameworks, simulation types, and simulation algorithms.

## Programmatically retrieving information about simulation tools via the REST API

An API is available for retrieving information about the simulation tools in the BioSimulators collection. Please see the [documentation](https://api.biosimulators.org/) for the REST API for more information.

An interactive tutorial for using the API is available from [Binder](https://tutorial.biosimulators.org/).
# Source model repositories

In addition to projects contributed by investigators, BioSimulations also contains projects created from several primary repositories including [BiGG](http://bigg.ucsd.edu/), [BioModels](http://www.ebi.ac.uk/biomodels/), [ModelDB](http://modeldb.science/), the [Physiome Model Repository](https://models.physiomeproject.org/), the [Resource Balance Analysis Model Repository](https://github.com/SysBioInra/Bacterial-RBA-models), [RuleHub](https://github.com/RuleWorld/RuleHub), and the [VCell Published Models Database](https://vcell.org/vcell-published-models). Prior to incorporation into BioSimulations, BioSimulations extensively standardizes, quality-controls, and debugs the models and simulations in these resources; generates missing simulations and reports of their results; adds additional visualizations for simulation results; standardizes and fills in missing metadata; and packages these resources into standardized archives.

<div class="logos">
<div class="logos-row">
    <a href="http://bigg.ucsd.edu/" rel="noopener" target="_blank" title="BiGG">
    <img class="zoom" src="/assets/images/about/partners/bigg.png" />
    </a>

    <a
    href="http://www.ebi.ac.uk/biomodels/"
    rel="noopener" target="_blank"
    title="BioModels"
    >
    <img
        class="zoom"
        src="/assets/images/about/partners/biomodels.svg"
    />
    </a>

    <!--
    <a
    href="https://cellcollective.org/"
    rel="noopener" target="_blank"
    title="Cell Collective"
    >
    <img
        class="zoom"
        src="/assets/images/about/partners/cell-collective.png"
    />
    </a>    

    <a href="http://ginsim.org/" rel="noopener" target="_blank" title="GINsim">
    <img class="zoom" src="/assets/images/about/partners/ginsim.svg" />
    </a>

    <a
    href="https://jjj.biochem.sun.ac.za/"
    rel="noopener" target="_blank"
    title="JWS Online"
    >
    <img class="zoom" src="/assets/images/about/partners/jws.svg" />
    </a>
    -->

    <a
    href="http://modeldb.science/"
    rel="noopener" target="_blank"
    title="ModelDB"
    >
    <img class="zoom" src="/assets/images/about/partners/modeldb.svg" />
    </a>

    <a
    href="https://models.physiomeproject.org/"
    rel="noopener" target="_blank"
    title="Physiome Model Repository"
    >
    <img class="zoom" src="/assets/images/about/partners/physiome.svg" />
    </a>

    <a
    href="https://rba.inrae.fr/models.html"
    rel="noopener" target="_blank"
    title="RBA"
    >
    <img class="zoom" src="/assets/images/about/partners/rba.png" />
    </a>

    <a
    href="https://vcell.org/vcell-published-models"
    rel="noopener" target="_blank"
    title="VCell Published Models Database"
    >
    <img class="zoom" src="/assets/images/about/partners/vcell.svg" />
    </a>
</div>
</div>
