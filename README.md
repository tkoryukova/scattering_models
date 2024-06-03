### Brief description:

Implementation of the Extream Scattering Event (ESE) models in `julia` and `Python` following papers:

- [_Fiedler+1994_](https://ui.adsabs.harvard.edu/abs/1994ApJ...430..581F/abstract) Fiedler, R., Dennison, B., Johnston, K. J., Waltman, E. B., & Simon, R. S. (1994). **A summary of extreme scattering events and a descriptive model.** _The Astrophysical Journal_, 430, 581-594.
- [_Clegg+1998_](https://iopscience.iop.org/article/10.1086/305344/meta) Clegg, A. W., Fey, A. L., & Lazio, T. J. W. (1998). **The Gaussian plasma lens in astrophysics: refraction.** _The Astrophysical Journal_, 496(1), 253.
- [_Dong+2018_](https://academic.oup.com/mnras/article/481/2/2685/5090977) Dong, L., Petropoulou, M., & Giannios, D. (2018). **Extreme scattering events from axisymmetric plasma lenses.** _Monthly Notices of the Royal Astronomical Society_, 481(2), 2685-2693.

### Installation:

```bash
$ git clone https://github.com/tkoryukova/scattering_models.git
$ cd scattering_models
```

### Structure:

- `julia` module with scattering models from _Clegg+1998_ and _Dong+2018_ (`scattering_models.jl`)
- `julia` module with root finding utils (`root_finding_utils.jl`)
- `Python` script with _Fiedler+1994_ model (`Fiedler1994_model.py`)
- three scripts that plots the pictures for three scattering models considered (`fiedler1994_model.py`, `plot_clegg1998.jl`, `plot_dong2018.jl`)

### Example light curves:

![alt text](https://github.com/tkoryukova/scattering_models/blob/master/Fiedler1994.png?raw=true)
![alt text](https://github.com/tkoryukova/scattering_models/blob/master/Clegg1998.png?raw=true)
![alt text](https://github.com/tkoryukova/scattering_models/blob/master/Dong2018.png?raw=true)
