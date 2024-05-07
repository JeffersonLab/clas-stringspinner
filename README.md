# clas-stringspinner
[StringSpinner](https://gitlab.com/albikerbizi/stringspinner) wrapper for CLAS

## Dependencies

- [Pythia 8](https://pythia.org/)
- [fmt](https://fmt.dev/)
- [meson](https://mesonbuild.com/)
- [ninja](https://ninja-build.org/)

StringSpinner is included as a `git` submodule.

## Building

Build using `meson`. For example, let
- `$source_dir` be the `clas-stringspinner` source code directory (contains this `README.md`)
- `$build_dir` be the build directory
- `$install_dir` be the installation destination, which must be an absolute path

Build and install with:
```bash
meson setup $build_dir $source_dir --prefix $install_dir
meson install -C $build_dir
```

## Build options

| Option                          | Description                       |
| ---                             | ---                               |
| `stringspinner:install_example` | Install the StringSpinner example |

For more details, run:
```bash
meson configure $build_dir
```

To set a build option named `opt` to value `val`, run:
```bash
meson configure $build_dir -Dopt=val
```
