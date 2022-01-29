# Developer notes

## Technology stack

- [Single Page App](https://en.wikipedia.org/wiki/Single-page_application): The `cffinit` app will be a SPA. So app feels like a native app and no server-side code needs to run.
- [Figma](https://www.figma.com/): A vector graphics and prototyping editor used to developed the wireframes and interaction designs.
- [npm CLI](https://docs.npmjs.com/cli/v7): Package manager command line interface shipped with NodeJS.
- [TypeScript](https://www.typescriptlang.org/): Typed JavaScript language used for lowering maintenance cost.
- [Vue.js v3](https://v3.vuejs.org/): A frontend JS framework for building user interfaces.
- [Vue.js Composition API](https://v3.vuejs.org/guide/composition-api-introduction.html): Is style of writing UI components to group logical concerns like state management.
- [Quasar](https://quasar.dev/): A UI framework, a group of UI components, that follow [Material design guidelines](https://material.io/design) to make a coherent/pleasing user interface.
- [Quasar CLI](https://quasar.dev/quasar-cli): Build and test infrastructure. Combines well with Vue.js, Quasar and jest.
- [GitHub pages](https://pages.github.com/): Hosting of static HTML files. The build app (in `docs` folder) is deployed on it.
- [Husky](https://typicode.github.io/husky/#/): Automaticly runs checks before pushing changes to GitHub.
- [Jest](https://jestjs.io/): Testing framework to run unit tests and perform test assertions.
- [highlight.js](https://highlightjs.org/): To syntax highlight the YAML formatted file.
- [ESLint](https://eslint.org/): To get constistent code style and prevent errors the industry standard linter ESLint is used.

The notes about how we came to this technology stack, design and personas can be found in [project-docs/](project-docs/) folder.

## Clone the repository

```shell
# clone this repository
git clone https://github.com/citation-file-format/cff-initializer-javascript
# change directory
cd cffinit
```

## install dependencies

The command below will install dependencies

```shell
npm clean-install
```

## start the development server

```shell
npm run dev
```

Use a browser to navigate to [localhost:8080](http://localhost:8080/) to see the website.

## build the application

The command below will build the application and save the output in `docs/` folder.

```shell
npm run build
```

## linting the code

```shell
npm run lint
```

try to automatically fix linting issues with

```shell
npm run lint -- --fix
```

To run linting on commit, you can install a git commit hook with

```shell
npx husky install
```

## Tests

We use Jest for unit tests. To run unit tests (`test/jest/__tests__/**/*.jest.spec.ts`)

You can run the test with

```shell
npm run test:unit:ci
```

You can also use the Majestic web interface to run the unit tests in your browser.

```shell
npm run test:unit:ui
```

## Making a release

This section describes how to make a release in 2 parts:

1. preparation
1. making a release on GitHub

### (1/2) Preparation

1. Verify that the information in `CITATION.cff` is correct
2. Generate an updated version of `.zenodo.json` if needed using `cffconvert`
3. Make sure the version field in `package.json` is correct 
4. By running `npm run lint` make sure the linter does not complain
5. Run the unit tests with `npm run test:unit:ci`
6. Make sure that github.io page is up to date
7. Check whether the [Publish](https://github.com/citation-file-format/cff-initializer-javascript/actions/workflows/publish.yml) workflow worked recently and it was successful

### (2/2) GitHub

Make a [release on GitHub](https://github.com/citation-file-format/cff-initializer-javascript/releases/new).
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1404735.svg)](https://doi.org/10.5281/zenodo.1404735)
 
# cffinit: a web form to initialize CITATION.cff files

- Check out the **live version** [here](https://citation-file-format.github.io/cff-initializer-javascript/).
- For the rationale behind CITATION.cff files, read [the blog](https://www.software.ac.uk/blog/2017-12-12-standard-format-citation-files).
- For the Citation File Format specification, go [here](https://github.com/citation-file-format/citation-file-format) (latest) or [here](https://doi.org/10.5281/zenodo.1003149) (stable).
- For the Citation File Format home page, go [here](https://citation-file-format.github.io).
## Landing
![](00-landing.png)

## Start
![](01-start.png)

## Authors
![](02-authors.png)

## Identifiers
![](03-identifiers.png)

## Related resources
![](04-related-resources.png)

## Abstract
![](05-abstract.png)

## Keywords
![](06-keywords.png)

## License
![](07-license.png)

## Version specific
![](08-version-specific.png)

## Finish
![](09-finish.png)

## User stories

- Maria works for a funder en wants to ensure the funding helps as many people as possible. For a paper which depends on software for its results, that means you have to be able to accurately cite the software used to produce the results produced in the paper. She wants to make it obligatory for newly funded research to make sure their code is not only open but also citable, but she feels that she can only recommend this once she's confident that the website will give the users a good experience
- Shelly is a publisher at AGU. She wants to maximize the value of a given paper. Papers that are **reproducible** because they accurately cite what software was used for a paper, are more valuable.
- Peter is a developer who contributed to an existing project and now wants to also be **recognized** as an author of the software. The software already has a CITATION.cff file. He wants to be able to **load** the CITATION.cff information somehow, and then just add his own details and be done.
- Bouwe is a developer for a software package that has **many contributors**. His project does not yet have a CITATION.cff file. He wants to **automatically populate** the CITATION.cff file with the information about whoever is listed as a **contributor on his repository**.
- Stefano is a developer for a software package that has **many contributors**. His project does not yet have a CITATION.cff file but has been published on Zenodo. He wants to **automatically populate** the CITATION.cff file with the **information associated with the DOI**.
- Jan is a Phd student, who has developed some code that he wants people to cite. He has stumbled upon a link in the GitHub docs that point to cffinit. He wants to get this done in 5 minutes. Jan wants an easy way of getting the CITATION.cff file from cffinit to his repo. It should be clear to him that the minimal amount of work that is required of him is in fact a small amount of work, and he can come back later to make it better.
- Julia is a developer whose software was cited, but incorrectly because the citers didn't have access to **accurate** metadata pertaining to the software. She would like to have way to supply future citers with unambiguous data regarding e.g. the name of the software, the list of authors, etc
- Yoda is a Professor. He doesn't develop much software anymore himself, but more junior members of his research group do. He wants his group to be credited for the work they put in. He wants to send out an email to everybody in his group explaining why they need to make their software citable, and wants to include a link to cffinit to help people get started.
- Jezza is a researcher who wants to be able to accurately cite a piece of software that they **did not develop themselves**


## out of scope

- citing software developed by commercial companies


## Requirements

### milestones

1. Replicate the green cffinit v1 site's behavior but nicer
    - relevant user stories: Jan  
3. Add editing options based on existing CITATION.cff in a given github repo that users supply
    - relevant user stories: Peter 
5. Two way editing form -> text, text -> form
6. Add bringing in information from GitHub API, and offer functionality to merge the two information sources
    - relevant user stories: Bouwe
8. Add bringing in information from Zenodo API, and offer functionality to merge the two information sources
    - relevant user stories: Stefano
10. Add support for preferred-citation
11. Add support for references





# Pull request details

## List of related issues or pull requests

Refs: #ISSUE_NUMBER


## Describe the changes made in this pull request

<!-- include screenshots if that helps the review -->


## Instructions to review the pull request

<!--

```shell
cd $(mktemp -d --tmpdir cffinit-pr.XXXXXX)
git clone https://github.com/citation-file-format/cff-initializer-javascript .
git checkout <this branch>
npm clean-install
npm run dev
# go to localhost:8080, see if the app works correctly
npm run lint
npm run test:unit:ci
```

-->
