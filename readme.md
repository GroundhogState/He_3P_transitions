# Spectroscopy with He* BECs

This is a MATLAB suite for analysis of metrological data from the ANU He* lab.

**Purpose:**
This software uses particle detection instances and diagnostic metadata to determine portions the electronic spectrum of 4He\* from experimental data obtained in the ANU He\* BEC laboratory.

**Design philosophy:**
Favours readability and detailed diagnostics at some marginal expense.
Caching is implemented at important breakpoints after expensive operations

**Status:** This project is in working order, though functionality has not been tested exhaustively. Data can be made available upon reasonable request.


### Prerequisites

This code is written for MATLAB >= R2018a.
You'll need the [Core BEC analysis](https://github.com/brycehenson/Core_BEC_Analysis) library, which should be cloned into the same directory for ease of use.
NOTE: The txy_import function will not execute properly with some Unix distributions. It requires a creation time for data files, which is stored by default in Windows. This software is fully functional in Windows. 

### Installing

Clone into the directory of your choice, along with the required libraries.

### Getting started

The core program is main_scripts, which calls a series of functions to import, condition, and process data stored in our archive as comma-delimited .txt files. One day, this may become publically available.

## Contributing

Contributions are welcome, if you happen to have any use for this code, by the usual fork/edit/pull request.

## Authors

**Jacob Ross** [groundhogstate](https://github.com/groundhogstate)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* A substantial amount of prior work was done by [**Bryce Henson**](https://github.com/brycehenson).
* Thanks to [**Kieran Thomas**](https:/github.com/KF-Thomas) for assistance during data collection.

### Built With
* PurpleBooth's [readme template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* Our [Core BEC analysis](https://github.com/brycehenson/Core_BEC_Analysis), including several contributions to MATLAB file exchange


## TODO

* Clean up code that isn't in the dependency tree
* Document processing pipeline
* Add project logo
* Unit tests
* Upload data to AARNet