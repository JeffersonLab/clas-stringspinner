# clas-stringspinner
[StringSpinner](https://gitlab.com/albikerbizi/stringspinner) wrapper for CLAS

## Dependencies

- [Pythia 8](https://pythia.org/)
- [fmt](https://fmt.dev/)
- [meson](https://mesonbuild.com/)
- [ninja](https://ninja-build.org/)

StringSpinner is included as a `git` submodule.

## Building

Build using `meson`, _e.g._:

```bash
meson setup build /path/to/clas-stringspinner --prefix /path/to/installation/prefix
meson install -C build
```
