# Crystal Mountain

A script for regression testing of command line executables that output text (E.g. compilers).

All it does is run small command line files that end in `.case.sh`, manage their outputs and detect changes. It's a
very simple Python program but has proven to be useful. It was first written to test a compiler built for a university course.

## Getting started

 1. Copy, submodule or otherwise get `crystalmountain.py` to be callable from a directory in your project
 2. In such directory, create a `cases` subdirectory
 3. Populate `cases` with `*.case.sh` files that run the program you want to test

And you are al set. The suite is operational. You can run `crystalmountain.py -g` to generate initial outputs. Which 
will be written to corresponding `*.spec` files.

Note that cases can also be in nested subdirectories

### Examples
Assuming the above setup:

 * `crystalmountain.py`: run the test suite. This will record if there are any changes in the outupts
 * `crystalmountain.py cases/my/test.case1.sh cases/my/test.case2.sh`: run those two particular test cases
 * `crystalmountain.py -g`: generate new outputs, will not override existing ones
 * `crystalmountain.py -s`: show the diffs from the last run
 * `crystalmountain.py -l -f`: list all the test cases that are failing


### Full Usage
```text
usage: crystalmountain.py [-h] [-p | -c | -s | -r | -g | -l] [-f | -n | -a] [CASE ...]

Run a regression test suite for liblsqecc

positional arguments:
  CASE                Specify one or more particular test cases

optional arguments:
  -h, --help          show this help message and exit
  -p, --displayonly   display the results of the execution without writing anything
  -c, --copyover      Instead of running, copy over the results from the previous run
  -s, --showdiffs     Show diffs after failed tests
  -r, --removeoutput  Removes outputs for this test
  -g, --generate      Generates tests if they don't exist yet, otherwise do nothing
  -l, --list          List available tests
  -f, --failing       Only consider cases that failed in the latest run
  -n, --nospec        Only consider cases that don't have a spec defined
  -a, --passing       Only consider cases that passed in the last run

```

### Requirements
The requirements are pretty minimal:

 * A shell to run commands (E.g. bash or zsh)
 * A Python >= 3.7 installation
 * A `diff` command for rendering changes (Optional)

Only tested on Linux systems, but probably works on MAC WSL/GitBash.

## License

Details of the license in the top of `crystalmountain.py`
