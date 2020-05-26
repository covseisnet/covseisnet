__Roadmap__
===


## General

- __Authors__: how do we decide to mention the codes authors? As the code is going to evolve, it is going to be interesting to have specific authors per functions, with associated email adress.

- __Docstrings__: according to book [_The Hitchhiker's Guide to Python: Best Practices for Development_](https://docs.python-guide.org/writing/documentation/), we should develop the code according the [PEP8 standards](https://www.python.org/dev/peps/pep-0008/), which is what the Black soft is doing to your code if these standards are not met. Within the same book, [sphinx](https://www.sphinx-doc.org/en/master/) is presented as the most popular documentation tool that allow to generate a documentation webpage from the codes' docstrings. Is is based on the reStructured text formatting used in the docstring. The reStructured text format is quite flexible, and several docstring standards exist. In order to be as simple as possible, most of the code so far uses the [numpy docstring style](https://numpydoc.readthedocs.io/en/latest/format.html#). Should we agree on this?

- __Documentation__: based on the above comments, the documentation is generated from the code docstring. After installing `sphinx`, we should also install the following sub-packages:
    - `sphinx_bootstrap_theme` for webpage style
    - `sphinx_execute_code` for example execution on build

    The main advantage of using `sphinx` is for allowing inter-sphinx mapping for the documentation (referencing a function based on the alread-existing documentation of other libraries like numpy, etc. that are also coded in sphinx).

- __Setup__: for installing the package, we have to declare what is required. I made a non-exhaustive list in the freshly created [setup.py](setup.py) (a copy from this [setup template](https://github.com/pypa/sampleproject/blob/master/setup.py)).

- verbose mode ? `logger`, `tqdm`
- top-level function from examples?
    - obspy example coherence
    - 1-y long script on PdF
- everything is a numpy array (np.save/z, np.load)

## Data management

- Remove methods for strange formats (`h5read`, `matread`). We want the users to whether feed covariance with an `obspy.Stream` object stream or an `covseisnet.ArrayStream` object.

- for calculating the covariance, if we get a `obspy.Stream` object, we have to turn it into a `covseisnet.ArrayStream`, involing a `synchronize` routine (method of the `data.py` module).

- the `time` propertry returns only 1 time vector from 1 of the traces, knowing that they should be exactly on the same sampling and length. If not, this property should trigger an error/warning.

- the `homegenize` will become `synchronize` for synchronizing streams on same sampling and length automatically (with many options).

- stop using `savitzky_golay` for smoothing. Move to `scipy.signal.smooth`.

- stop using robust for claculating the `mad`. Move to `scipy.signal.mad`.

- remove `"method"` ``kwargs, when smooth=1, it is onebit spectral whitening.

- removed plot method.


## Covariance

- move fft to covariance
- remark: fft really depends on parameters
- francis if you want to try differrnt fft up to you
- `spectral_width` et `entropy` sont dans `coherence`

### Todo
- Cyril parle de covariance non-hermitienne. Il faut définir les propriétés d'une telle matrice (comment dealer avec les vp non hermitiennes, comment on calcule la coherence dans ce cas).
- Spectral lines
- `rank` peut être une liste de ranks.
- inclure iso et formulate synthetic according to it

## Correlation

### Todo
- Sigma en secondes


## Net is the new antenna

- FDSN
- Exemples: montrer les correlation sur Kamchatka triées avec net.epicentral_dist_from.


## Beamforming (planewave)

-

## Location ("beam" sphérique)

### Calcul ttimes

```
cov, net = np.load('cov.npz')
# cov.select()
# net.select()
cov = cov.eigenvectors(rank=0)
corr = cn.correlation.calculate(cov)

model = np.load('model.npz')
# dictionary with 'lon', 'lat', 'z'

"ttimes.npz"
ttimes = np.load(tties.npz)
ttimes[station1]
ttimes[station2]


# Define grid
likelihood = covseisnet.GridSearch(model, net, ttimes='/path/to/ttimes.npz', build=False)
likelihood.calculate(corr)
```

- On prévoit de coder la sélection + tard (Francis!)
