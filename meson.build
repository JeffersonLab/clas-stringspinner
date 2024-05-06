project(
  'clas-string-spinner',
  'cpp', 'fortran',
  meson_version: '>=1.2',
  default_options: {
    'cpp_std':   'c++17',
    'buildtype': 'release',
    'libdir':    'lib',
  },
  version: '0.0.0',
)

# build stringspinner
run_command(
  meson.project_source_root() / 'subprojects/stringspinner/build.sh',
  get_option('clean_stringspinner') ? 'clean' : '',
  check: true,
  capture: false,
)
