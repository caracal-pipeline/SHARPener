# SHARPener

This is a set of tools that have been developed in preparation of the [SHARP survey](
https://www.astron.nl/astronomy-group/apertif/science-projects/sharp-search-hi-absorption-apertif/sharp). 

The main function of `sharpener` is to identify the position of all continuum sources in a continuum image and extract
a spectrum from each line of sight of these sources. 

The spectra are then plotted. 

`sharpener` can be run using a `.yml` parameter file (link to default) as `python sharpener.py <parameter_file.yml>`, or through a `IPython`
notebook. 

The following tutorials can guide you through the different capabilities of `sharpener`.

### Installation

**Requisites**
- SHARPener makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`). 
- The parameter file is in `yaml` format, hence `pyaml`, and `json` packages should be installed

**Insallation instructions**
- Clone this repository. From terminal type:

```
git clone https://github.com/Fil8/SHARPener.git
```

- add `sharpener` directory to `PYTHONPATH` (for permanent use save it in your `.cshrc_profile`, or `.bash_profile`, respectively)

```
setenv PYTHONPATH $path_to_sharpener:${PYTHONPATH}

export PYTHONPATH=$PYTHONPATH:path_to_sharpener
```
