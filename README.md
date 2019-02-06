# SHARPener


This is a set of tools that have been developed in preparation of the [SHARP survey](
https://www.astron.nl/astronomy-group/apertif/science-projects/sharp-search-hi-absorption-apertif/sharp). 

The main function of `sharpener` is to identify the position of all continuum sources in a continuum image and extract
a spectrum from each line of sight of these sources. 

The spectra are then plotted. 

`sharpener` can be run automatically using a `.yml` [parameter file](https://github.com/Fil8/SHARPener/wiki/Parameter-file) as `python sharpipeline.py <path_to_parameter_file.yml>`, or through a `IPython`
[notebook](https://github.com/Fil8/SHARPener/blob/master/tutorials/T2_automated_run.ipynb). 

The following [tutorials](https://github.com/Fil8/SHARPener/tree/master/tutorials) can guide you through the different capabilities of `sharpener`.

***

### Installation

**Requisites**
- SHARPener makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`) and addition to `mpdaf` and `pypdf2`. 
- The parameter file is in `yaml` format, hence `pyaml`, and `json` packages should be installed.
- Tutorials make use of `tabulate` and `glob` for fancy outputs.

**Insallation instructions**
- Clone this repository. From terminal type:

```
git clone https://github.com/Fil8/SHARPener.git
```

- add `sharpener` directory to `PYTHONPATH` (for permanent use save it in your `.cshrc_profile`, or `.bash_profile`, respectively)

```
setenv PYTHONPATH path_to_sharpener:${PYTHONPATH}

export PYTHONPATH=$PYTHONPATH:path_to_sharpener
```

- change path at `line16` of `sharpener.py` as follows: `sys.path.append('/path-to-sharpener/sharp_modules/')` 
 
 ***
 <p>&copy <sub> Filippo M. Maccagni 2018-2019 </sub></p>
