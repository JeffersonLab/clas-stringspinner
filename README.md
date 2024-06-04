# clas-stringspinner
[StringSpinner](https://gitlab.com/albikerbizi/stringspinner) wrapper for CLAS

## Dependencies

- [Pythia 8](https://pythia.org/)
- [fmt](https://fmt.dev/)
- [meson](https://mesonbuild.com/)
- [ninja](https://ninja-build.org/)
- [StringSpinner](https://gitlab.com/albikerbizi/stringspinner) is already included as a `git` submodule

## Building

Build using `meson`. For example, let
- `$source_dir` be this `clas-stringspinner` source code directory (contains this `README.md`)
- `$build_dir` be the build directory (will be created)
- `$install_dir` be the installation destination, which must be an absolute path (will be created)

1. Setup the build directory:
```bash
meson setup $build_dir $source_dir --prefix $install_dir
```
2. Set build options (see **Build Options** section below); the defaults should be sufficient for most users, so this step is **optional**.
3. Build and install:
```bash
meson install -C $build_dir
```
4. Test:
```bash
meson test -C $build_dir
```

### Build options

In addition to the built-in options, the following table shows the project options:

| Option                          | Description                                                                                            |
| ---                             | ---                                                                                                    |
| `stringspinner:install_example` | Install the [StringSpinner example](https://gitlab.com/albikerbizi/stringspinner/-/blob/master/dis.cc) |

For more details and current option values, run
```bash
meson configure $build_dir
```
To set a build option named `opt` to value `val`, run:
```bash
meson configure $build_dir -D opt=val
```

## Running StringSpinner

The installed binary `stringspinner` is the main executable Monte Carlo generator; assuming you have not changed certain build options, run
```bash
$install_dir/bin/stringspinner --help
```
for the usage guide.

> [!TIP]
> For easier usage, add the installed `$install_dir/bin` directory to your `$PATH`, or make a symbolic link pointing to the `stringspinner` executable.
