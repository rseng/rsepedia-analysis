# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

--------------------------------------------------------------------------------------------------------------------------------

## [Unreleased]

### Added
- paper.bib - Citation details for Zenodo archive of v1.0.0
- README.md - Citation details for JOSS paper


### Changed
- README.md - Updated Development Status section

### Fixed
- paper.md - A few typos and formatting issues

### 

## [v1.0.0] - 2020-08-29 - First Official Release

### Added
- .travis.yml - Build matrix entries and corresponding script code for testing clang compiler
- docs/User_Manual.pdf - Copy of user manual to the docs directory for easy linking and viewing on the web
- makefile - Detection of versioned g++ compilers
- makefile - Detection of clang compiler and setting the compiler flags
- msvc - Directory with Microsoft Visual studio solution and project files for building Excimontec on Windows
- paper/paper.md - short paper describing Excimontec for publication in JOSS
- paper/paper.bib - Bibtex reference information for the JOSS paper
- README.md - Recommended reading, citation information, and acknowledgement sections
- README.md - Information about clang compiler support
- User_Manual.pdf - User manual PDF to the root directory including examples for each type of simulation test, a detailed description of the input parameters and the phyiscal models implemented in the code, and detailed installation instructions
- user_manual - Directory with source files for the user manual including all of the simulation output data used for the examples

### Changed
- .gitignore - Ignore statements for Microsoft Visual Studio files to not ignore the solution and project files but still ignore the build directories
- .travis.yml - Updated copyright statement years
- .travis.yml - Updating testing config by removing testing of GCC v4.7 and v4.8 and added testing of GCC v9
- KMC_Lattice - KMC_Lattice submodule to v2.1.0
- LICENSE - Updated copyright statement years
- main.cpp - Version string to v1.0.0
- main.cpp - Updated copyright statement years
- makefile - Updated googletest directory to the one located within the KMC_Lattice submodule
- makefile - Updated copyright statement years
- makefile - Updated how the underlying compiler is detected now using 'mpicxx -show'
- parameters_default.txt - Version string to v1.0.0
- parameters_default.txt - Parameter section headings
- README.md - Replaced detailed installation and build instructions with link to new user manual
- README.md - Reorganized sections to be more useful for new users
- README.md - Updated information about Travis CI testing configuration
- README.md - Updated copyright statement years

### Removed
- .travis.yml - Build matrix entries for testing GCC 4.x
- .travis.yml - Build matrix entries for mpich with GCC 6 and 8
- .gitmodules - The googletest submodule entry
- .travis.yml - Coveralls exlcude statement for the googletest directory
- .travis.yml - sudo statement because it is no longer used by Travis CI
- googletest - Duplicate googletest submodule with the intent to use the googletest submodule already within the KMC_Lattice submodule

### Fixed
- .travis.yml - Changed Linux build matrix entries to use xenial or bionic instead of trusty to fix code coverage reporting error

## [v1.0.0-rc.3]- 2019-04-01 - Density of States Integration Bugfix

### Added
- OSC_Sim - Steady_DOS_sampling_count and Steady_DOOS_sampling_counter to keep track of how many times the DOS and DOOS are sampled during the simulation
- test.cpp (SteadyTransportTests) - Tests to check the integral of the DOS and DOOS with and without including Coulomb interactions

### Changed
- main.cpp - Version string to v1.0.0-rc.3
- OSC_Sim (updateSteadyData) - To increment the new DOS and DOOS sampling counters and removed calculation of the DOS because the DOS is does not change during the simulation when not including Coulomb interactions
- OSC_Sim (getSteadyDOS) - Function is not longer const
- parameters_default.txt - Version string to v1.0.0-rc.3

### Fixed
- main.cpp - Output of density of states and density of occupied states data to have the correct units
- OSC_Sim (getSteadyDOOS) - To use the new sampling counter to avoid errors averaging over multiple samplings
- OSC_Sim (getSteadyDOOS_Coulomb) - To use the new sampling counter to avoid errors averaging over multiple samplings
- OSC_Sim (getSteadyDOS) - To calculate the DOS using only one sample at call time because the DOS (without Coulomb) does not change during the simulation
- OSC_Sim (getSteadyDOS_Coulomb) - To use the new sampling counter to avoid errors averaging over multiple samplings

## [v1.0.0-rc.2]- 2019-02-09 - Steady State Charge Transport Test Update

### Added
- .gitignore - Numerous ignore statements to ignore files generated during build and test operations and from Microsoft Visual Studio
- CHANGELOG.md - Notes about all changes in this release
- docs - KMC_Lattice documentation
- docs - Documentation of new functions in the OSC_Sim class
- Exciton - Project name (Excimontec) to header guards
- main.cpp - Output of equilibration energy and transport energy both with and without including the Coulomb potential when running a steady transport test
- main.cpp - Output the current density when running a steady transport test
- main.cpp - Output the density of states and density of occupied states data calculated with and without Coulomb potential adjustments
- OSC_Sim - Project name (Excimontec) to header guards
- OSC_Sim - Public function and member variable documentation
- OSC_Sim (getSiteEnergy) - Checking of the input coordinates validity and generating error if invalid
- OSC_Sim (getSiteType) - Checking of the input coordinates validity and generating error if invalid
- OSC_Sim - include statements for all used components to better document class dependencies
- OSC_Sim (getSteadyCurrentDensity) - New function that returns the average current density from the steady state transport test
- OSC_Sim (getSteadyEquilibrationEnergy_Coulomb) - New function to return the equilibration energy calculated including the Coulomb potential
- OSC_Sim (getSteadyTransportEnergy_Coulomb) - New function to return the transport energy calculated including the Coulomb potential
- OSC_Sim (Steady_equilibration_energy_sum_Coulomb) - New private member variable for calculating the equilibration energy including the Coulomb potential
- OSC_Sim (Transport_energy_weighted_sum_Coulomb) - New private member variable for storing the data used to calculate the transport energy including the Coulomb potential
- OSC_Sim (Steady_hops_per_DOS_sample) - New private member variable for calculating the DOS during the steady transport test
- OSC_Sim (Steady_hops_per_DOOS_sample) - New private member variable for calculating the DOOS during the steady transport test
- OSC_Sim (DOS_bin_size) - New private member variable for calculating the DOOS and DOS during the steady transport test
- OSC_Sim (executePolaronHop) - Calculation of the transport energy including the Coulomb potential
- OSC_Sim (exportEnergies) - New overloaded function allowing the user to output the absolute site energies for electrons or holes
- OSC_Sim (generateSteadyPolarons) - Simpler creation of polarons on random sites when energetic disorder is disabled
- OSC_Sim (generateSteadyPolarons) - Status output about how many polarons are created in the lattice for the steady transport test
- OSC_Sim (generateSteadyPolarons) - Allow creation of polarons on acceptor sites
- OSC_Sim (getSteadyDOOS) - New function for getting density of occupied states data
- OSC_Sim (getSteadyDOOS_Coulomb) - New function for getting density of occupied states data
- OSC_Sim (getSteadyDOS) - New function for getting density of states data
- OSC_Sim (getSteadyDOS_Coulomb) - New function for getting density of states data
- OSC_Sim (getSteadyTransportEnergy) - Function now allows calculation when the sum of weights is negative
- OSC_Sim (getSteadyTransportEnergy_Coulomb) - Function now allows calculation when the sum of weights is negative
- OSC_Sim (outputStatus) - Status output for the steady state charge transport test
- OSC_Sim (updateSteadyData) - Calculation of the equilibration energy including the Coulomb potential
- OSC_Sim (updateSteadyData) - Periodic sampling of the DOOS and DOS
- OSC_Sim (updateSteadyDOS) - New function that uses the input energy of a given site to update the input density of states data
- Parameters - Project name (Excimontec) to header guards
- Parameters - Public function and member variable documentation
- Parameters (Enable_steady_data_output) - New member variable that is set to true by default but could be used in the future to disable DOS data output
- Polaron - Project name (Excimontec) to header guards
- Polaron - Public function and member variable documentation
- README.md - Build instructions link for Windows users
- README.md - Information about new DOS and DOOS data files generated during the steady transport test
- README.md - Examples about what the custom site energies import feature can be used for.
- test.cpp - Added a simple command line status message at the beginning of all test cases
- test.cpp (SteadyTransportTests) - Tests checking the output of the transport and equilibration energies when the steady transport test has not been run
- test.cpp (SteadyTransportTests) - Tests comparing the energies calculated with and without Coulomb interactions and relative positions of the equilibration and transport energies
- test.cpp - (SteadyTransportTests) - Test to check that phase restriction disabling increases the number of available sites for creating the initial polarons in donor-acceptor blends
- test.cpp - (SteadyTransportTests) - Tests to check the peak position of the DOS and DOOS data from the very low field test
- test.cpp (SteadyTransportTests) - Tests to check the transport energy calculated during a medium electric field test condition
- test.cpp (SteadyTransportTests) - Test to check the relative position of the transport energy and the donor HOMO during the medium electric field test
- test.cpp (SteadyTransportTests) - Test to check the absolute position of the transport energy
- test.cpp (SteadyTransportTests) - Test for transport in a random donor-acceptor blend and check for the relative change in transport energy position relative to the neat donor architecture
- test.cpp (SteadyTransportTests) - Test checking that the simulation works when phase restriction in enabled with a random blend
- test.cpp (SteadyTransportTests) - Test checking the magnitude of the current density
- test.cpp (ToFTests) - Test checking the hole extraction map output

### Changed
- KMC_Lattice - KMC_Lattice submodule to v2.1.0-beta.1
- Many files - Copyright statement years to 2017-2019
- docs - Updated docs using Doxygen v1.8.15
- Doxyfile - Project version number to v1.0.0
- Doxyfile - Settings so that KMC_Lattice is included and markdown files are no longer be included in the documentation
- Exciton - Nested derived Exciton event classes into the Exciton class
- Exciton (constructor) - Must now specify the spin state upon construction
- Exciton - Revised documentation for public functions and member variables
- main.cpp - Version string to v1.0.0-rc.2 in preparation for next release
- OSC_Sim - Functions to use new object event class nesting format
- OSC_Sim (Site_OSC) - Store the site type internally as a char instead of short to save memory
- OSC_Sim (generateExciton) - Moved definition of default tag value to the header file function declaration statement
- OSC_Sim (generateExciton) - Implemented the new Exciton constructor where one must specify the spin state
- OSC_Sim - Nested derived Site class (Site_OSC) into the OSC_Sim class as a private class
- OSC_Sim (calculateDOSCorrelation) - Scope of the two functions from public to private
- OSC_Sim (executePolaronHop) - Refactored function and using temporary local variables to make code more readable
- OSC_Sim (executePolaronHop) - Calculation of the transport energy to absolute value that includes the HOMO energy and accounts for donor and acceptor site occupation
- OSC_Sim (executePolaronHop) - Calculation of the transport energy is only performed after the equilibration phase is complete
- OSC_Sim (executePolaronHop) - Calculation of the transport energy is is done using the displacement as the weights instead of the velocity
- OSC_Sim (getChargeExtractionMap) - Refactored code replacing usage of stringstream with addition of substrings
- OSC_Sim (updateSteadyData) - Calculation of equilibration energy to absolute value that includes the HOMO energy and accounts for donor and acceptor site occupation
- Parameters (checkParameters) - Lower the limit for the smallest internal potential that one can during a steady transport simulation to allow very low field simulations
- parameters_default.txt - Default morphology file format to not include the compression specifier suffix 
- Polaron - Nested derived Polaron event classes into the Polaron class
- README.md - Replaced version badges with text links
- README.md - Updated description of steady transport test feature
- README.md - Updated current release status info for KMC_Lattice to v2.1
- test.cpp - All tests to route command line output to a test_log.txt file instead of cluttering the command line making it easier to see the test results
- test.cpp (EnergiesImportTests) - Tests of the new exportEnergies function checking the absolute value of the exported electron and hole site energies
- test.cpp (ExcitonDynamicsTests) - Increased the range of the transient to get a more accurate assessment of the equilibrium energy position
- test.cpp (LoggingTests) - Adjusted parameters to promote additional mechanisms to occur including RISC and exciton recombination
- test.cpp (SteadyTransportTests) - Test of the energy values to compare to absolute energy values including the HOMO energy
- test.cpp (SteadyTransportTests) - Very low field test to make it more accurate and a little bit faster by decreasing the internal potential, lattice size, and Coulomb cutoff radius
- test.cpp (SteadyTransportTests) - Medium field test by reducing the number of iterations to make the test faster
- test.cpp (SteadyTransportTests) - Adjusted parameters of the no disorder mobility test to decrease the test time

### Removed
- docs - Markdown files from the generated documentation
- main.cpp - Output of the Fermi energy during the steady state charge transport test
- OSC_Sim (Steady_Fermi_energy) - private member variable that is no longer used
- OSC_Sim (getSteadyFermiEnergy) - Fermi energy is no longer calculated and DOOS and DOS data is output instead
- test.cpp (SteadyTransportTests) - Tests of the Fermi energy

### Fixed
- CHANGELOG.md - Several spelling mistakes in previous release sections
- main.cpp - Spelling mistake in the results output of steady transport test
- main.cpp - Label error in the results output of the time-of-flight charge transport test, current should have been current density
- main.cpp - Bug bug where status output and logfile reset was only being performed when checking the error status
- OSC_Sim (calculateCoulomb) - Specified input parameter namespaces to avoid Doxygen confusion between the header declaration and source file definition
- OSC_Sim (createExciton) - Specified input parameter namespaces to avoid Doxygen confusion between the header declaration and source file definition
- OSC_Sim (calculatePolaronEvents) - Spelling mistake in error message
- OSC_Sim (reassignSiteEnergies) - Spelling mistake in error message
- OSC_Sim (executePolaronHop) - Transport energy calculation to correctly account for hops across the periodic boundary when calculating the displacement
- OSC_Sim (getSteadyTransportEnergy) - Spelling mistake in the function documentation
- OSC_Sim (getPolaronIt) - Bug that could occur when comparing a list iterator to an iterator from a different list
- OSC_Sim (getSteadyEquilibrationEnergy) - Bug that could cause rounding error due to integer division
- Parameters (checkParameters) - Spelling mistakes in error messages
- Parameters (importParameters) - Spelling mistakes in error messages
- parameters_default.txt - Incorrect units for the exciton hopping and annihilation rate prefactors
- test.cpp - Spelling mistakes in the test comments

## [v1.0.0-rc.1]- 2018-12-11 - Interfacial Energy Shift, Site Energy Import, and Steady State Transport Test Update

### Added
- CHANGELOG.md - New file detailing the changes for each release
- README.md - Link to new Changelog file
- README.md - Copyright statement
- slurm_script.sh - Copyright statement
- Parameters - New class to store all parameters used by the simulation and to contain functions for parsing the parameter file and checking for parameter validity
- Parameters - Default values to parameters used in main
- makefile - Commands to compile and link the Parameters class source files
- makefile - Missing dependency of all Excimontec objects on the KMC_Lattice lib
- test.cpp (ParameterTests) - New test for importing the parameters_default.txt file that checks the importParameters function
- test.cpp (ParameterTests) - New series of tests for how the program handles misspelled boolean values in the parameter file
- parameters_default.txt - Three new parameters (Enable_interfacial_energy_shift, Energy_shift_donor, Energy_shift_acceptor) for the interfacial energy shift model
- Parameters - The three new parameters for the interfacial energy shift feature
- Parameters (importParameters) - Code to read the interfacial energy shift parameters from the parameter file
- Parameters (checkParameters) - Validity checks for the interfacial energy shift parameters
- OSC_Sim (reassignSiteEnergies) - Code to implement the interfacial energy shift model
- test.cpp - The three new parameters for the interfacial energy shift model to the default parameters struct
- test.cpp (InterfacialEnergyShiftTests) - New test function checking the energy shift using a bilayer with and without energetic disorder
- test.cpp (ParameterTests) - Tests to check response of OSC_Sim to initialization with invalid interfacial energy shift parameters
- parameters_default.txt - Two new parameters (Enable_import_energies, Energies_import_filename) for importing site energies from a file
- Parameters - The two new parameters for the site energy import feature
- Parameters (importParameters) - Code to read the site energy import parameters from the parameter file
- Parameters (checkParameters) - Validity checks for the new site energies import parameters
- OSC_Sim (exportEnergies) - New function that creates a site energies file
- OSC_Sim (reassignSiteEnergies) - Code section that imports the energies from the specified file when the site energy import feature is enabled
- New site energies files that have improper format to check for handling of invalid site energies files during testing
- test.cpp - The two new site energies import parameters to the default parameters struct
- test.cpp - (ParameterTests) - Tests to check response of OSC_Sim to initialization with invalid site energies import parameter combinations
- test.cpp (EnergiesImportTests) - New test function checking the export and import of valid energies file and to check how the program handles energies file with improper format or missing data
- Feature allowing users to enable event logging for debugging purposes by adding -enable_logging as a command line argument after the parameter file name
- main.cpp (main) - Code to check command line arguments and enable logging if the -enable_logging argument is passed after the parameter file name
- OSC_Sim - New N_events_executed counter and getN_events_executed function to overwrite the one in the base Simulation class to separately keep track oh how many events have been executed
- OSC_Sim - New private member variable previous_event_type to store the type of event that was just executed and a getPreviousEventType function to retrieve the string for testing purposes
- OSC_Sim (calculateNextEvent) - Code that increments the new N_events_executed and store the previous event type
- test.cpp (LoggingTests) - New test function to check that all events being executed are accurately logged and test the getN_events_executed as well
- test.cpp (ChargeDynamicsTests) - New simple test for an exciton dissociation and charge recombination dynamics simulation
- test.cpp (ExcitonDiffusionTests) - New test to check that increased energetic disorder reduces the diffusion length, which also tests activated exciton hopping functionality
- test.cpp (ExcitonDiffusionTests) - New tests checking the exciton recombination event counters and a test for exciton diffusion in an exponential DOS
- test.cpp (IQETests) - New test with a weakly donating/accepting bilayer to test its impact on the charge separation yield and allow polarons to hop to opposing phase
- test.cpp (ObjectCreationTests) - New test function to test the createExciton, createElectron, and createHole functions of the OSC_Sim class with both valid and invalid input coordinates and site types and to test some of the event counters and the checkFinished function
- test.cpp (ParameterTests) - Several previously untested invalid parameter combinations
- test.cpp (ToFTests) - New hole ToF test using the Marcus hopping model
- test.cpp (ToFTests) - New tests to check that increased energetic disorder reduces the mobility, which also tests activated polaron hopping functionality
- test.cpp (ToFTests) - New test that checks the behavior of the ToF_placement_energy feature
- test.cpp (ToFTests) - New tests that get and analyze the transient energy relaxation data and test various aspects of this behavior
- README - New information about the interfacial energy shift model, the site energy import capability, and the steady transport test
- makefile - New PGI compiler flag enabling output of compiler warnings
- parameters_default.txt - New parameters for the steady transport test: Enable_steady_transport_test, Steady_carrier_density, and N_equilibration_events
- main.cpp (main) - Command line output and results file output for the steady transport test
- OSC_Sim - New public functions for the steady transport test: getSteadyEquilibrationEnergy, getSteadyFermiEnergy, getSteadyMobility, and getSteadyTransportEnergy
- OSC_Sim- New private member variables to store data needed to calculate the final steady transport test results: Steady_equilibration_energy_sum, Steady_equilibration_time, Transport_energy_weighted_sum, and Transport_energy_sum_of_weights
- OSC_Sim - New generateSteadyPolarons function for creating initial polarons a the beginning of the steady transport test
- OSC_Sim - New updateSteadyData function to updating the steady transport test data variables each event iteration
- OSC_Sim (init) - Call to the new generateSteadyPolarons function when the steady transport test is enabled
- OSC_Sim (calculatePolaronEvents) - Code to determine when polarons are trying to cross the z periodic boundary and adjust the potential energy change accordingly
- OSC_Sim (calculatePolaronEvents) - Code to disable calculation of polaron extraction events when the steady transport test is enabled
- OSC_Sim (checkFinished) - Code to check for the test termination conditions for the steady transport test
- OSC_Sim (executeNextEvent) - Call to the new updateSteadyData function before executing each event when the steady transport test is enabled
- OSC_Sim (executePolaronHop) - Code to record the transport energy data when the steady transport test is enabled
- OSC_Sim (updateTransientData) - Code to cast the return of the distance function as an int to prevent compiler warnings
- Parameters - New parameters for the steady transport test: Enable_steady_transport_test, Steady_carrier_density, and N_equilibration_events
- Parameters (checkParameters) - Call to the base class checkParameters function
- Parameters (checkParameters) - New checks for the steady transport test parameter combinations
- Parameters (importParameters) - Code to import new steady transport test parameters from the parameter file
- Parameters (importParameters) - Code to check the status of the input ifstream and throw an exception if there is a problem
- test.cpp - Default parameter values for the new parameters for the steady transport test: Enable_steady_transport_test, Steady_carrier_density, and N_equilibration_events
- test.cpp (ParameterTests) - New tests checking the validity of the parameters when enabling the steady transport test
- test.cpp - New test function SteadyTransportTests that check the output of the steady transport test
- test.cpp (ParameterTests) - New tests for attempting to import parameters using an uninitialized or close ifstream
- OSC_Sim - New overloaded createExciton function that does not specify the creation coordinates and generates them randomly
- OSC_sim (createExciton(Coords)) - Code to check whether the input coordinates are occupied and producing an error if they are
- test.cpp (ChargeDynamicsTests) - New tests to check the transient data
- test.cpp (ExcitonDynamicsTests) - New tests to check exciton relaxation into the tail of a Gaussian DOS, triplet dissociation, and triplet annihilation mechanisms
- test.cpp (IQETests) - New tests for changes in the number of geminate recombination events under different test conditions
- test.cpp (ObjectCreationTests) - New tests checking attempts to create an exciton in a fully occupied lattice
- test.cpp (ParameterTests) - New tests for invalid parameters in the Parameters_Lattice and Parameters_Simulations classes
- test.cpp (SteadyTransportTests) - New tests for attempts to create more initial polarons that there are donor sites
- test.cpp (ToFTests) - New tests for attempts to create more initial polarons that there are donor sites
- test.cpp (ToFTests) - New tests comparing the relaxed mobility and relaxed occupation energy from the transient data against the expected steady state mobility and steady state equilibration energy
- CONTRIBUTING.md - New file with detailed instructions for how new people can contribute to the project
- README.md - Link to the CONTRIBUTING.md file

### Changed
- main.cpp (main) - Refactored code to use the new Parameters class
- OSC_Sim - Refactored all functions to use to use the new Parameters class member
- test.cpp - Refactored all tests to use the new Parameters class
- test.cpp - Renamed ParameterCheckTests to ParameterTests
- OSC_Sim (getSiteEnergy) - Updated function access from private to public so that it can be used during testing
- OSC_Sim (getSiteType) - Updated function access from private to public so that it can be used during testing
- parameters_default.txt - Formatting of interfacial energy shift model parameters
- .travis.yml - Copyright to show correct years, 2017-2018
- .travis.yml - Switched test coverage compiler configuration from GCC 5 to GCC 4.7 to make the test faster
- LICENSE - Copyright to show correct years, 2017-2018
- makefile - Copyright to show correct years, 2017-2018
- OSC_Sim (calculateExcitonCreationCoords) - Refactored code to be more robust and be guaranteed to find an appropriate empty site for creation if one exists
- KMC_Lattice - Submodule to latest version that fixes the chooseNextEvent and removeObject bugs that occurred when an object does not have a valid event
- OSC_Sim (createElectron, createExciton, CreateHole) - Code to catch out_of_range exceptions from the lattice class as the method for checking for invalid input coordinates
- test.cpp (EnergiesImportTests) - Test to use a bilayer so that it tests assignment of site energies to both donor and acceptor type sites
- test.cpp (ExcitonDiffusionTests) - Updated test by increasing N_tests to gather more statistics for checking the numerical accuracy of the lifetime and reducing N_tests when simply checking for relative changes in the diffusion length
- test.cpp (ExcitonDynamicsTests) - Updated test to use reduced Dynamics_transient_end and check how the program handles cutting off the tail end of the transient
- test.cpp (ExcitonDynamicsTests) - Updated test by reducing N_tests to shorten the test time
- Doxyfile - Updated the release number to v1.0.0-rc.1
- KMC_Lattice - Updated to latest version that has additional features needed by the steady transport test
- parameters_default.txt - Updated the version number to v1.0.0-rc.1
- main.cpp (main) - Updated version string to v1.0.0-rc.1
- main.cpp (main) - Updated Parameters object usage to use new format from the updated KMC_Lattice submodule
- Site_OSC - Storage of site energies using pointers to now directly storing a float
- OSC_Sim - Updated all appropriate functions to use float data type for site energies instead of double to save memory
- Parameters - Update base class from the struct used by the old KMC_Lattice submodule version to the new Parameters_Simulation class in the latest KMC_Lattice submodule version
- Parameters (importParameters) - Updated function to use format from new Parameters_Simulation base class in the latest KMC_Lattice version
- test.cpp - Updated functions to use format from new Parameters_Simulation base class in the latest KMC_Lattice version
- test.cpp (InterfacialEnergyShiftTests) - Updated site energies to use float data type instead of double
- test.cpp (ExcitonDynamicsTests) - Increased the lifetime tests tolerance to reduce likelihood of test failure
- main.cpp (main) - Updated results file output to use new getN_singlet_excitons_dissociation and getN_triplet_excitons_dissociated functions instead of the old getN_excitons_dissociated function
- OSC_Sim (calculateExcitonCreationCoords) - Changed function name to calculateRandomExcitonCreationCoords
- OSC_Sim (getN_excitons_dissociated) - Replace the function with getN_singlet_excitons_dissociated and getN_triplet_excitons_dissocated allowing the user to get separate stats
- OSC_Sim (N_excitons_dissociated) - Replace the counter variable with N_singlet_excitons_dissociated and N_triplet_excitons_dissociated to keep separate stats
- OSC_Sim (executeExcitonDissociation) - Changed event counter by counting the singlet and triplet dissociation events separated using the new counter variables
- test.cpp - Reduced default test parameters ToF_pnts_per_decade and Dynamics_pnts_per_decade down to 10 from 20 to have a larger time step to average over
- test.cpp (SteadyTransportTests) - Adjusted N_equilibration_events and N_tests to reduce test time
- README.md - Updated status to note that all features are now implemented and have major testing
- docs - Updated Doxygen documentation

### Removed
- main.cpp (Parameters_main) - the Parameters_main struct because these parameters are now stored on the Parameters class
- main.cpp (importParameters) - the importParameters function because this function is now contained in the Parameters class
- OSC_Sim (Parameters_OPV) - the Parameters_OPV struct because these parameters are now stored in the Parameters class
- OSC_Sim - Many of the parameter variables (replaced by a Parameters class member)
- OSC_Sim (checkParameters) - the checkParameters function because this function is now contained in the Parameters class
- OSC_Sim (createImportedMorphology) - Several unused local variables
- OSC_Sim (getEvents) - The function had no immediate purpose
- Parameters (checkParameters) - Several parameter checks that are now performed by the base Parameters_Simulation class
test.cpp (ParameterTests) - Several tests that are now handled by the KMC_Lattice submodule

### Fixed
- OSC_Sim (calculatePolaronEvents) - Update code so an error is not generated if no valid events are calculated because sometimes polarons can be trapped on sites where there are no valid events
- Parameters (checkParameters) - Bug where the function was not indicating an error when unable to read in the bool value for the Enable_recalc option
- Parameters (checkParameters) - Bug where the check for enabling of both imported energies and interfacial energy shift was implemented wrong
- OSC_Sim (getChargeExtractionMap) - Bug where the function as skipping output of the 0,0 coordinate data

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.5] - 2018-11-01 - Morphology Import and BKL Algorithm Update

### Added
- README.md - More detailed instructions for building and testing Excimontec
- CODEOWNERS.md - New file that designates which users have ownership over which parts of the codebase
- CODE_OF_CONDUCT.md - New file that specifies standards for developer and user interactions
- Excimontec namespace for all source code
- Usage of new KMC_Lattice namespace throughout code (replaces Utils namespace)
- main.cpp (main) - Usage of the new calculateTransitTimeHist function from the OSC_Sim class
- makefile - Commands to the clean target, so that make clean is also run for the KMC_Lattice submodule
- OSC_Sim - morphology files are imported so they are now opened and checked by the OSC_Sim class for easier testing
- Exciton.h - New calculateRateConstant function to the derived Event classes
- OSC_Sim (createImportedMorphology) - Code to allow import of uncompressed morphologies from Ising_OPV v3.2 and v4.0
- OSC_Sim (createImportedMorphology) - Code that parses and checks the morphology file version and handles attempts to import invalid morphology files
- test.cpp (ToFTests) - Electron ToF simulation tests
- test.cpp (ExcitonDynamicsTests) - Triplet exciton lifetime dynamics test
- test.cpp (ExcitonDiffusionTests) - Triplet exciton diffusion test
- test.cpp (IQETests) - Several basic IQE tests for field, temperature, delocalization, and recombination rate dependence of the charge separation yield in a bilayer morphology
- test.cpp (IQETests) - IQE testing for high illumination intensity to check whether higher order exciton-exciton and exciton-polaron annihilations and bimolecular recombination events increase as expected
- test.cpp (ParameterCheckTests) - Many new tests checking for how the OSC_Sim class handles initialization with invalid parameters
- Several new valid and invalid test morphology files from Ising_OPV to be used during morphology file import tests
- test.cpp (MorphologyImportTests) - New test function MorphologyImportTests and added numerous tests of importing valid and invalid sample morphology files to check for correct behavior

### Changed
- docs - Updated Doxygen API documentation
- README.md - Section and text formatting 
- README.md - Name for the transit time histogram file from "ToF_transit_time_dist.txt" to "ToF_transit_time_hist.txt"
- KMC_Lattice - Submodule to latest development branch version that has the BKL algorithm implemented, the Version class, and other bugfixes
- Code formatting for all source files using MS Visual Studio style
- Copyright in files to show correct years, 2017-2018
- main.cpp (main) - Name of output file for transit time data from "ToF_transit_time_dist.txt" to "ToF_transit_time_hist.txt" to specify that it is histogram data
- main.cpp (main) - Usage of calculateTransitTimeDist function to use new calculateTransitTimeHist that was renamed in the OSC_Sim class
- main.cpp (importParameters) - Usage of the new str2bool function from the KMC_Lattice namespace
- OSC_Sim (calculateExcitonEvents, calculatePolaronEvents) - Refactored code to use new BKL algorithm from the KMC_Lattice library
- OSC_Sim - Replaced morphology file ifstream with a morphology filename string
- OSC_Sim (calculateTransitTimeDist) - Function to output vector of pairs to match the format of the rest of the software package
- OSC_Sim (calculateTransitTimeDist) - Renamed function to calculateTransitTimeHist to be more descriptive about the data that it produces
- OSC_Sim (generateDynamicsExcitons) - Reduced frequency of command line status output for dynamics simulations
- OSC_Sim (generateToFPolarons) - Reduced frequency of command line status output for ToF simulations
- OSC_Sim (generateDynamicsExcitons) - Dynamics simulation command line output to print info when starting the first transient cycle
- OSC_Sim (executeExcitonHop, executePolaronHop) - Command line error output to give more details about the event when there is an error
- .travis.yml - Method for running the test executable
- test.cpp (ToFTests) - Refactored test to use updated calculateTransitTimeHist function
- test.cpp (MorphologyImportTests) - Test morphology file names to include their path relative to the root directory so that tests can be run from the root directory
- test.cpp (ExcitonDynamicsTests) - Reduced N_tests to shorten the test calculation time

### Removed
- Usage of Utils namespace throughout code (replaced by KMC_Lattice namespace)
- main.cpp (Parameters_main) - Morphology filename member variable
- main.cpp (main) - Morphology file ifstream
- main.cpp (main) - Opening and checking morphology files for import 
- main.cpp (importParameters) - Usage of the importBooleanParam function that was removed from the KMC_Lattice library
- Exciton.h - removed calculateExecutionTime overloaded definition from all derived Event classes
- Polaron.h - removed calculateExecutionTime overloaded definition from all derived Event classes
- OSC_Sim - Stopped compile warning about undefined pragma lines by commenting the out the pragma lines

### Fixed
- OSC_Sim (checkParameters) - A few parameter check bugs

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.4] - 2018-05-17 - Testing and Continuous Integration Update

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.3] - 2018-02-16

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.2] - 2018-02-08

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.1] - 2017-10-23
<!---
# Copyright (c) 2017-2020 Michael C. Heiber
# This source file is part of the Excimontec project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec
--->
# Excimontec

Kinetic Monte Carlo simulations are a powerful computational tool that have been used in concert with experiments and more detailed theoretical methods to understand and optimize organic semiconductor materials and devices. 
However, despite over 30 years of applying KMC tools to organic semiconductors, no widespread or standardized software tools have taken hold in the community. 
Instead, many research groups around the world have maintained private codebases of varying complexity, efficiency, and reliability. 
As a result, there have been large barriers to entry for new researchers and a lot of repeated effort throughout the community that would be much better off applied to pushing the capabilities of the technique and further refining the physical models. 

Excimontec represents an honest effort to bring the community together around a well-tested, optimized, reliable, and accessible open-source tool for performing KMC simulations of organic electronic devices. 
The software is being developed in modern C++ and is optimized for efficient execution on high performance computing clusters using MPI. 
This software package uses object-oriented design and extends the [KMC_Lattice](https://github.com/MikeHeiber/KMC_Lattice) framework. 
If you would like to contribute to the development of this project, please see the [contributing instructions](CONTRIBUTING.md).

Major releases and other significant developments will be announced on the Excimontec: General News mailing list. If you are interested in keeping up to date with the latest developments, please subscribe at the following link:
[Subscribe Here](http://eepurl.com/dis9AT)

#### Major Features:
- Adjustable periodic boundary conditions in all three directions allow users to perform 1D, 2D, or 3D simulations.
- Choose between several film architectures, including a neat film, bilayer film, or random blend film.
- Import bulk heterojunction morphologies generated by [Ising_OPV](https://github.com/MikeHeiber/Ising_OPV) v3.2 and v4.
- Donor and acceptor materials can take on an uncorrelated Gaussian density of states, a correlated Gaussian density of states with different correlation functions, or an uncorrelated exponential density of states model.
- Site energies at the donor-acceptor interface or in mixed regions can be modified using an interfacial energy shift model to generate an energy cascade.
- Custom site energies can also be imported from a text file to allow more exotic DOS distributions or implement the electrostatic potential due to fixed ionic dopants.
- Dynamics test simulations can be performed to generate exciton and charge carrier density transients that can be used to model exciton dissociation, charge carrier separation, and charge carrier recombination kinetics.
- Time-of-flight charge transport simulations of electrons or holes can be performed on neat, random blend, or bulk heterojunction blend films.
- Steady state charge transport simulations of holes can be performed on neat, random blend, or bulk heterojunction blends films with 3D periodic boundaries to estimate quasi-equilibrium properties including the steady state mobility, transport energy, and equilibration energy under different test conditions. Output of density of states and density of occupied states data can also be used to calculate the Fermi energy.  The Fermi energy and transport energy can then be used to calculate the Seebeck coefficient to investigate thermoelectric properties.
- Exciton diffusion simulations can be performed on any film architecture.
- Internal quantum efficiency simulations can be performed on bilayer, random blend, or bulk heterojunction blend films.
- Simulate complex exciton dynamics with events for intersystem crossing between singlet and triplet states as well as exciton-exciton and exciton-polaron annihilation events.
- Choose between Miller-Abrahams or Marcus models for polaron hopping.
- Charge carrier delocalization can be modeled with a spherical Gaussian delocalization model.
- Choose between several KMC algorithms (first reaction method, selective recalculation method, or full recalculation method).

For more details and examples of each type of simulation test, please see the [User Manual](https://mikeheiber.github.io/Excimontec/User_Manual.pdf)

## How to try Excimontec?

#### Building and Testing the Executable

This software tool uses [Message Passing Interface (MPI)](https://computing.llnl.gov/tutorials/mpi/) to utilize parallel computing power. 
As a result, using Excimontec requires that an MPI library is pre-installed on your system, and the final Excimontec executable must be built on your specific system. 
We cannot provide pre-built binaries for your system. 
Contact your HPC admin to determine the protocols for building MPI applications on your HPC system. 
In many cases, the HPC system will already be configured for you, and Excimontec comes with a default makefile that can be used with the [GCC compiler](https://gcc.gnu.org/), the [PGI compiler](https://www.pgroup.com/), or the [clang compiler](https://clang.llvm.org/).
If your system uses another compiler, you will need to edit the makefile and define your own compiler options.

If you wish, you can also install MPI on your own personal workstation and then build Excimontec there as well. 
For development and preliminary simulation tests, sometimes it is more efficient to run on your own workstation instead of an HPC system. 

For detailed installation instructions, please see the [User Manual](https://mikeheiber.github.io/Excimontec/User_Manual.pdf#page=13).

Please report any build or testing errors in the [Issues](https://github.com/MikeHeiber/Excimontec/issues) section. 

#### Usage

In most cases, your HPC system will use a job scheduler to manage the computing workload. 
For performing Excimontec simulations, it is recommended to submit batch jobs where you will request the resources needed to perform the simulation. 
An example batch script for the [SLURM](https://slurm.schedmd.com/) job scheduling system is provided with this package (slurm_script.sh). 
Similar batch scripts can also be written for other job schedulers.

Regardless of the job scheduler, the program execution command is essentially the same. 
Excimontec.exe takes one required input argument, which is the filename of the input parameter file. 
An example parameter file is provided with this package (parameters_default.txt).

For example, within the batch script, to create a simulation that runs on 10 processors, the execution command is

```mpiexec -n 10 Excimontec.exe parameters_default.txt```

In this example, the parameters_default.txt file that is located in the current working directory is loaded into the Excimontec program to determine which simulation to run.

#### Output

Excimontec will create a number of different output files depending which test is chosen in the parameter file:
- results#.txt -- This text file will contain the results for each processor where the # will be replaced by the processor ID.
- analysis_summary.txt -- When MPI is enabled, this text file will contain average final results from all of the processors.
- dynamics_average_transients.txt -- When performing a dynamics test, calculated exciton, electron, and hole transients will be output to this file.
- ToF_average_transients.txt -- When performing a time-of-flight charge transport test, calculated current transients, mobility relaxation transients, and energy relaxation transients will be output to this file.
- ToF_transit_time_hist.txt -- When performing a time-of-flight charge transport test, the resulting polaron transit time probability histogram will be output to this file.
- ToF_results.txt -- When performing a time-of-flight charge transport test, the resulting quantitative results are put into this parsable delimited results file.
- Charge_extraction_map#.txt -- When performing a time-of-flight or IQE test, the x-y locations where charges are extracted from the lattice are saved into this map file.
- DOS_correlation_data#.txt -- When a correlated density of states model is enabled, data showing the statistical correlation of site energies vs. distance is output into this file.
- DOS_data.txt - When performing a steady state charge transport test, the density of states distribution is calculated and output to this file.
- DOS_Coulomb_data.txt - When performing a steady state charge transport test, the density of states distribution, where the site energies include the Coulomb potential from the polarons, is calculated and output to this file.
- DOOS_data.txt - When performing a steady state charge transport test, the density of occupied states distribution is calculated and output to this file.
- DOOS_Coulomb_data.txt - When performing a steady state charge transport test, the density of occupied states distribution, where the site energies include the Coulomb potential from the polarons, is calculated and output to this file.

#### Data Analysis

For [Igor Pro](https://www.wavemetrics.com/) users, I am developing an open-source procedures package for loading, analyzing, and plotting data from Excimontec simulations called [Excimontec_Analysis](https://github.com/MikeHeiber/Excimontec_Analysis). 
This is a good starting point for managing the data generated by Excimontec, and the Igor Pro scripting environment provides a nice playground where users can perform more advanced data analysis as needed.

## Contact

If you would like some help in using or customizing the tool for your research, please contact me (heiber@mailaps.org) to discuss a collaboration. 
You can check out my research using this tool and other work on [Researchgate](https://www.researchgate.net/profile/Michael_Heiber).

Have a quick question or want to chat about Excimontec?  Join the discussion on Gitter: [![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg?style=for-the-badge)
](https://gitter.im/Excimontec)

## Citing This Work

[M. C. Heiber, J. Open Source Software **5**, 2307 (2020).](https://doi.org/10.21105/joss.02307) [[ResearchGate]](https://www.researchgate.net/publication/344141654_Excimontec_v10_An_Open-Source_Software_Tool_for_Kinetic_Monte_Carlo_Simulations_of_Organic_Electronic_Devices)

In addition, please also cite the DOI for the specific version that you used from [Zenodo.org](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%22595806%22&sort=-version&all_versions=True).

## Recommended Reading

Below are some recommended resources for starting to learn about KMC modeling of organic electronic devices:

[M. C. Heiber, A. Wagenpfahl, and C. Deibel, "Advances in Modeling the Physics of Disordered Organic Electronic Devices" In *Handbook of Organic Materials for Electronic and Photonic Devices*, Woodhead Publishing Series in Electronic and Optical Materials, edited by O. Ostroverkhova (Woodhead Publishing, 2019) 2nd Ed., Chap. 10.](https://doi.org/10.1016/B978-0-08-102284-9.00010-3) [[ResearchGate]](https://www.researchgate.net/publication/329625990_Advances_in_Modeling_the_Physics_of_Disordered_Organic_Electronic_Devices)

[C. Groves, "Simulating Charge Transport in Organic Semiconductors and Devices: A Review" Rep. Prog. Phys. **80**, 026502 (2017).](http://dx.doi.org/10.1088/1361-6633/80/2/026502)[[ResearchGate]](https://www.researchgate.net/publication/311743859_Simulating_charge_transport_in_organic_semiconductors_and_devices_A_review)

[S. D. Baranovskii, "Theoretical description of charge transport in disordered organic semiconductors" Phys. Status Solidi B. **3** 487 (2014).](https://doi.org/10.1002/pssb.201350339)

[C. Groves, "Developing Understanding of Organic Photovoltaic Devices: Kinetic Monte Carlo Models of Geminate and Non-geminate Recombination, Charge Transport and Charge Extraction" Energy Environ. Sci. **6**, 32020 (2013).](https://doi.org/10.1039/c3ee41621f)[[ResearchGate]](https://www.researchgate.net/publication/260303139_Developing_Understanding_of_Organic_Photovoltaic_Devices_Kinetic_Monte_Carlo_Models_of_Charge_Recombination_Transport_and_Extraction)

## Development Status

The current release, [Excimontec v1.0.0](https://github.com/MikeHeiber/Excimontec/releases), is built with [KMC_Lattice v2.1](https://github.com/MikeHeiber/KMC_Lattice/releases) and allows the user to perform several simulation tests relevant for OPV and OLED devices. 
Please report any bugs or submit feature requests in the [Issues](https://github.com/MikeHeiber/Excimontec/issues) section. 
Please see the [Changelog](CHANGELOG.md) for a detailed listing of previous and upcoming changes. 

#### Continuous Integration and Testing Status:

Excimontec is currently being tested on [Ubuntu](https://www.ubuntu.com/) (v16.04 and v18.04) with the [GCC compiler](https://gcc.gnu.org/) (v5, v6, v7, v8, and v9) the [clang compiler](https://clang.llvm.org/) (v7) and on both [Open MPI](http://www.open-mpi.org/) (v1.10.2) and [MPICH](http://www.mpich.org/) (v3.2 and v3.3) using [Travis CI](https://travis-ci.com/).

| Branch | Status |
| :------: | ------ |
| Master | [![Build Status](https://img.shields.io/travis/MikeHeiber/Excimontec/master.svg?style=for-the-badge)](https://travis-ci.org/MikeHeiber/Excimontec) |
| Development | [![Build Status](https://img.shields.io/travis/MikeHeiber/Excimontec/development.svg?style=for-the-badge)](https://travis-ci.org/MikeHeiber/Excimontec) |

Code is being tested using [googletest](https://github.com/google/googletest) with test coverage assessment by [Coveralls](https://coveralls.io/).

| Branch | Status |
| :------: | ------ |
| Master | [![Coveralls github branch](https://img.shields.io/coveralls/github/MikeHeiber/Excimontec/master.svg?style=for-the-badge)](https://coveralls.io/github/MikeHeiber/Excimontec?branch=master) |
| Development | [![Coveralls github branch](https://img.shields.io/coveralls/github/MikeHeiber/Excimontec/development.svg?style=for-the-badge)](https://coveralls.io/github/MikeHeiber/Excimontec?branch=development) |

## For Software Developers

Public API documentation for the Excimontec package can be viewed [here](https://mikeheiber.github.io/Excimontec/).

## Acknowledgments
Thank you to Dr. Dean M. DeLongchamp at NIST for providing access to computing resources that have supported the development of v1.0. 
Development of v1.0 has been primarily supported by financial assistance award 70NANB14H012 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD).
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

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
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
reported by contacting the project team at heiber@mailaps.org. All
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

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contributing to Excimontec

We welcome new contributors to this project to help make Excimontec better and more effective for our scientific and educational users. Contributions can take on a variety of forms, such as bug reports, feature requests, code documentation updates, user tutorials, feature implementations, bugfixes, and more.

## Ground Rules

- To maintain an inclusive and welcoming environment, all community members must follow the [Code of Conduct](./CODE_OF_CONDUCT.md).
- To keep project development organized, all new features and bugfixes must be proposed and discussed using [issues](https://github.com/MikeHeiber/Excimontec/issues) before changes are implemented.
- C++ code should follow the [C++11 standard](https://en.wikipedia.org/wiki/C%2B%2B11) and should be written to be simple and largely self-documenting with self-descriptive names.

## Filing Issues

One of the easiest ways to contribute is by using the software and submitting bug reports and feature requests as issues in the [Issues section](https://github.com/MikeHeiber/Excimontec/issues). 
Before filing a new issue please check that your issue has not already been submitted. 
If it has, 'like' the issue and add any additional feedback you have to the issue thread in a comment.
If your issue has not already been filed, create a new issue and provide a detailed description of the problem or new idea/suggestion. 
The code owners will review new issues, provide feedback, and decide if and when to implement the changes needed to resolve the issue.

## Working on an Issue

1. Before you set out to tackle an issue, discuss your implementation plans with the code owners on the issue page.
    
    Especially for code additions needed to add new features, it is important to discuss your implementation ideas with the code owners to ensure that the software design is relatively well thought out and fits with the rest of the codebase.
    
2. Create fork of the Excimontec repository on your account.

    Click the "Fork" button near the top of any Excimontec Github webpage.

3. Connect your fork to the upstream repository

    ```
    git remote add upstream https://github.com/MikeHeiber/Excimontec.git
    ```

4. Clone your forked Excimontec repository to your local machine.

    ```
    git clone --recurse-submodules https://github.com/YOUR_USERNAME/Excimontec
    ```

5. Create a new feature branch from the development branch.

    ```
    git checkout development
    git branch feature-new-thing
    git checkout feature-new-thing
    ```

6. Implement your changes keeping the scope of changes limited to the specified issue.

7. Document your changes by adding to the [Changelog](./CHANGELOG.md).

8. Update your new branch with any upstream changes and resolve any conflicts.

    ```
    git fetch upstream
    git checkout development
    git merge upstream/development
    git checkout feature-new-thing
    git rebase development
    ```

9. Make sure your changes do not break any of the existing functionality by running the tests locally.

    ```
    make test
    ./test/Excimontec_tests.exe
    ```
    
10. Push your local changes to your forked Excimontec repository.

    ```
    git push origin feature-new-thing
    ```

11. Submit a pull request

    See below for details about how to submit a pull request.

## Code Testing and Documentation Standards

Unit tests using the [googletest framework](https://github.com/abseil/googletest/blob/master/googletest/docs/primer.md) must be created for all new/modified functions and added to the test/test.cpp or test/test_mpi.cpp file as appropriate. Tests should include comments that describe what functionality each test code section is testing for.

All new/modified classes and functions must be documented using the [Doxygen style](https://www.stack.nl/~dimitri/doxygen/manual/docblocks.html).

## Submitting A Pull Request

The pull request must reference the issue that is being resolved.
Merging of feature branches into the development branch will be done using a squash merge operation to keep a tidy Git history.
The squash merge will be done by a code owner, so you must provide a single overall commit message that details all changes in the pull request body text for the code owner the use with the merge.
The pull request version of the code must pass all [TravisCI tests](https://travis-ci.org/MikeHeiber/Excimontec) and must not reduce code testing coverage on [Coveralls](https://coveralls.io/github/MikeHeiber/Excimontec).
Code will be formatted during the code review stage to the Microsoft Visual Studio 2017 standard format.
Updated documentation files will also be generated by Doxygen during the code review stage as needed.
Code owners may add new commits to your pull request branch in order to help you fix any problems.

---
title: 'Excimontec v1.0: An Open-Source Software Tool for Kinetic Monte Carlo Simulations of Organic Electronic Devices'
tags:
  - kinetic Monte Carlo
  - C++
  - organic photovoltaics
  - organic semiconductors
  - exciton diffusion
  - charge recombination
  - charge transport
  
authors:
 - name: Michael C. Heiber
   orcid: 0000-0002-1567-5663
   affiliation: "1"
affiliations:
 - name: Center for Hierarchical Materials Design (CHiMaD), Northwestern University, Evanston, Illinois 60208, USA
   index: 1
date: 4 May 2020
bibliography: paper.bib
---

# Summary

For over three decades, kinetic Monte Carlo (KMC) simulations have been a powerful computational tool to help understand and optimize organic semiconductor devices, especially photovoltaics, light-emitting diodes, transistors, and thermoelectrics [@baranovskii2014pssb; @groves2017rpp; @heiber2019chapter; @zuo2019aem].
KMC simulations use mechanistic models for how excitons and polarons are created, migrate through, and are then eventually removed from the semiconductor layer of a device and can capture the complex interactions between performance and spatial structure that is often not possible using continuum drift-diffusion models. 
This can then be used to probe a wide variety of phenomena in organic electronic devices, including exciton diffusion and quenching, charge transport, and charge recombination at the full device scale while retaining details regarding nanoscale inhomogeneities. 
Despite the clear utility of the method, no widespread or standardized software tools have taken hold in the community. 
Instead, many research groups around the world have maintained private codebases of varying complexity, efficiency, and reliability. 
As a result, there have been large barriers to entry for new researchers and a lot of repeated effort throughout the community that would have been much better off applied to pushing the capabilities of the technique and further refining the physical models.

``Excimontec`` is designed to be a well-tested, optimized, reliable, and accessible open-source tool for performing KMC simulations of organic electronic devices. 
v1.0 has a particular focus on organic photovoltaic device modeling and can utilize complex bulk heterojunction morphologies generated using the ``Ising_OPV`` tool [@heiber2018joss].
v1.0 comes with five different simulation tests: exciton diffusion, time-of-flight charge transport, internal quantum efficiency, dynamics, and steady state charge transport. 
See the user manual for a more in-depth description each simulation test, including some examples of what each test can be used for. 
The software has been developed in modern C++ and is optimized for efficient execution on high performance computing clusters using MPI. 
This software package uses object-oriented design and extends the ``KMC_Lattice`` framework [@heiber2019joss].
The code includes rigorous unit and validation testing with ``googletest``, continuous integration testing with ``TravisCI``, and API documentation generated using ``Doxygen``. 
The source code for ``Excimontec v1.0`` is archived with Zenodo [@heiber2020excimontec1.0.0].


# Acknowledgments

This work was developed under the financial assistance award 70NANB14H012 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD).  Thank you to Dr. Dean M. DeLongchamp for providing access to NIST's Raritan computing cluster, which was helpful with software development and testing.

# References
