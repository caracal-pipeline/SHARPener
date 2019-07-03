# SHARPener


This is a set of tools that have been developed in preparation of the HI absorption surveys of the SKA precursors and pathfinders.

`sharpener` identifies the position of all continuum sources in a continuum image and extracts a spectrum against all their lines of sight. The spectra are then plotted. Multiple options can be provided by the user, such has hanning smoothing and polynomial fitting and subtraction of the spectra, and plotting the spectra in different units (km/s, Hz, Jy beam$^{-1}$, $\tau$). Further information are given in the next sections.

`sharpener` can be run automatically using a `.yml` [parameter file](https://github.com/Fil8/SHARPener/wiki/Parameter-file) as  `run_sharpener -c <parameter_file.yml>`, or through a `IPython`
[notebook](https://github.com/Fil8/SHARPener/blob/master/tutorials/T2_automated_run.ipynb).

To generate config with `run_sharpener -gd <parameter_file.yml>` and for help `run_sharpener -h`.

The following [tutorials](https://github.com/Fil8/SHARPener/tree/master/tutorials) can guide you through the different capabilities of `sharpener`.

***

### Installation

**Requisites**
- `SHARPener` makes use of the most common `python` packages (e.g. `numpy`, `scipy`, `astropy`, `astroquery`, `prettytable`) and addition to `mpdaf` and `pypdf2`.
- The parameter file is in `yaml` format, hence `pyaml`, and `json` packages should be installed.
- Tutorials make use of `tabulate` and `glob` for fancy outputs.
- `SHARPener` source finder calls the function `imsad` of Miriad. Miriad must be installed for the source finder to work properly. No python wrapper is needed, since it is already included in `SHARPener`.

**Insallation instructions**

- This package is available on [pypi](https://pypi.org/project/sharpener/), allowing:

```
pip install sharpener
```

- Alternatively, clone this repository. From terminal type:

```
git clone https://github.com/Fil8/SHARPener.git
```

- Then change directory into SHARPener and install:

```
cd SHARPener && pip install .

```

 
 ***
 <p>&copy <sub> Filippo M. Maccagni 2018-2019 </sub></p>
