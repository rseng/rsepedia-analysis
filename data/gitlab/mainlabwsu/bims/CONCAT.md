## Project
BIMS (Breeding Information Management System) is a module that enables breeders to store, manage, archive and analyze their public or private breeding program data.


## Requirement
 - Drupal 7.x
 - Mainlab Chado Loader 7.x

## Version
1.0.0

##Download
The BIMS module can be download from GitLab:

https://gitlab.com/mainlabwsu/bims

## Installation
After downloading the module, extract it into your site's module directory 
(e.g. sites/all/modules) then follow the instructions below:

  1. Enable the module by using the Drupal administrative interface: Go to: Modules, check 'BIMS' and save or by using the 'drush' command:
     ```
     drush pm-enable bims
     ```
     This will create all BIMS related tables in public schema and populate the tables with default values. It will also create directories for BIMS in Drupal public file directory.

## Documentation
Please visit [BIMS home page](https://www.breedwithbims.org) to learn more about this module.

## Problems/Suggestions
BIMS is still under active development. For questions or bug report, please contact the developers at the Main Bioinformatics Lab by emailing to: dev@bioinfo.wsu.edu
